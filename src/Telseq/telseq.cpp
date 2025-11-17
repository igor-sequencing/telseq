
#include <iostream>
#include <istream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include <istream>
#include <sstream>
#include <set>
#include <map>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <future>

#include "telseq.h"
#include "Timer.h"
#include "math.h"
#include "Util.h"
#include "prettyprint.h"

// Use htslib instead of BamTools for URL streaming support
#include <htslib/sam.h>
#include <htslib/hts.h>

//
// Getopt
//

#define PROGRAM_BIN "telseq"
#define AUTHOR "Zhihao Ding"

static const char *TELSEQ_VERSION_MESSAGE =
"TelSeq Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR ".\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n" ;

static const char *TELSEQ_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n\n"
"Usage: " PROGRAM_BIN " [OPTION] <in.1.bam> <in.2.bam> <...> \n"
"Scan BAM and estimate telomere length. \n"
"   <in.bam>                 one or more BAM files to be analysed. File names can also be passed from a pipe, \n "
"                            with each row containing one BAM path.\n"
"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should has only one column, \n"
"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified, \n"
"                            <in.bam> are ignored.\n"
"   -i, --bam-index=STR      URL to BAM index file (.bai). Required when using -T > 1 with S3/HTTP URLs. \n"
"                            For local files, index is auto-detected. For URLs, you must specify the index URL.\n"
"   -o, --output_dir=STR     output file for results. Ignored when input is from stdin, in which case output will be stdout. \n"
"   -j, --jobs=INT           max number of BAM files to process in parallel. 0 = auto-detect cores. default = 0.\n"
"   -T, --threads-per-file=INT  number of threads for processing reads within each BAM file (requires BAM index). default = 1.\n"
"   -H                       remove header line, which is printed by default.\n"
"   -h                       print the header line only. The text can be used to attach to result files, useful\n"
"                            when the headers of the result files are suppressed. \n"
"   -m                       merge read groups by taking a weighted average across read groups of a sample, weighted by \n"
"                            the total number of reads in read group. Default is to output each readgroup separately.\n"
"   -u                       ignore read groups. Treat all reads in BAM as if they were from a same read group.\n"
"   -U                       process ONLY unmapped reads and exit. Requires BAM index. Useful for quick unmapped read analysis.\n"
"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
"\nTesting functions\n------------\n"
"   -r                       read length. default = 100\n"
"   -z                       use user specified pattern for searching [ATGC]*.\n"
"   -e, --exomebed=STR       specifiy exome regions in BED format. These regions will be excluded \n"
"   -w,                      consider BAMs in the speicfied bamlist as one single BAM. This is useful when \n"
"                            the initial alignemt is separated for some reason, such as one for mapped and one for ummapped reads. \n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static StringVector bamlist;
    static std::string outputfile = "";
    static std::string exomebedfile = "";
    static std::string bamindexurl = ""; // URL to BAM index file for S3/HTTP streaming
    static std::map< std::string, std::vector<range> > exomebed;
    static bool writerheader = true;
    static bool mergerg = false;
    static bool ignorerg = false;
    static bool unmappedonly = false; // process only unmapped reads and exit
    static bool onebam = false; // whether to consider all bams as one bam
    static int tel_k= ScanParameters::TEL_MOTIF_CUTOFF;
    static std::string unknown = "UNKNOWN";
    static std::string PATTERN;
    static std::string PATTERN_REV;
    static unsigned int parallel_files = 0; // max number of BAM files to process in parallel (0 = auto)
    static unsigned int threads_per_file = 1; // number of threads per BAM file for read processing

}

// Global mutex for thread-safe output
static std::mutex output_mutex;

// Helper function to check if a path is a URL
bool isURL(const std::string& path) {
    return (path.find("http://") == 0 ||
            path.find("https://") == 0 ||
            path.find("ftp://") == 0 ||
            path.find("s3://") == 0);
}

static const char* shortopts = "f:o:i:k:z:e:r:p:j:T:HhvmuwU";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bamlist",		optional_argument, NULL, 'f' },
    { "output-dir",		optional_argument, NULL, 'o' },
    { "bam-index",      optional_argument, NULL, 'i' },
    { "exomebed",		optional_argument, NULL, 'e' },
    { "jobs",           optional_argument, NULL, 'j' },
    { "threads-per-file", optional_argument, NULL, 'T' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// combine counts in two ScanResults objects
void add_results(ScanResults& x, const ScanResults& y){

    x.numTotal += y.numTotal;
    x.numMapped += y.numMapped;
    x.numDuplicates += y.numDuplicates;
    x.n_exreadsExcluded += y.n_exreadsExcluded;
    x.n_exreadsChrUnmatched += y.n_exreadsChrUnmatched;
    x.n_totalunfiltered += y.n_totalunfiltered;

    // Ensure x vector is large enough to hold all data from y
    std::size_t max_telcount_size = std::max(x.telcounts.size(), y.telcounts.size());
    if (x.telcounts.size() < max_telcount_size) {
        x.telcounts.resize(max_telcount_size, 0);
    }

    std::size_t max_gccount_size = std::max(x.gccounts.size(), y.gccounts.size());
    if (x.gccounts.size() < max_gccount_size) {
        x.gccounts.resize(max_gccount_size, 0);
    }

    // Now x is properly sized, iterate through y's actual size and add to x
    for (std::size_t j = 0; j < y.telcounts.size(); ++j){
        x.telcounts[j] += y.telcounts[j];
    }

    for (std::size_t k = 0; k < y.gccounts.size(); ++k){
        x.gccounts[k] += y.gccounts[k];
    }

}

// merge results in result list into one
void merge_results_by_readgroup(

    std::vector< std::map<std::string, ScanResults> >& resultlist){
    std::vector< std::map<std::string, ScanResults> > mergedresultslist;
    std::map<std::string, ScanResults> mergedresults;

    for(size_t i =0; i< resultlist.size(); i++){
        const auto& rmap = resultlist[i];
        for(std::map<std::string, ScanResults>::const_iterator it= rmap.begin();
                it != rmap.end(); ++it){

            const std::string& rg = it ->first;
            const ScanResults& result = it -> second;

            if(mergedresults.find(rg) == mergedresults.end()){
                mergedresults[rg] = result;
            }else{
                add_results(mergedresults[rg], result);
            }
        }
    }
    mergedresultslist.push_back(mergedresults);
    resultlist = mergedresultslist;

}

// Process a single chromosome from a BAM file (thread-safe function)
// This function is designed to be called in parallel for different chromosomes
std::map<std::string, ScanResults> processChromosome(
    const std::string& bampath,
    int tid,
    const std::string& chr_name,
    bool isExome,
    bool rggroups,
    const std::map<std::string, ScanResults>& resultTemplate) {

    // Create a clean copy with zero counts (only copy structure, not accumulated data)
    std::map<std::string, ScanResults> chromResults;
    for (const auto& pair : resultTemplate) {
        ScanResults cleanResult;
        cleanResult.sample = pair.second.sample;
        cleanResult.lib = pair.second.lib;
        cleanResult.bam = pair.second.bam;
        chromResults[pair.first] = cleanResult;
    }
    std::map<std::string, std::vector<range>::iterator> lastfound;
    std::vector<range>::iterator searchhint;

    // Open BAM file for this thread
    samFile* sam_fp = sam_open(bampath.c_str(), "r");
    if (sam_fp == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: could not open BAM file in chromosome thread: " << bampath << std::endl;
        return chromResults;
    }

    // Read header
    sam_hdr_t* header = sam_hdr_read(sam_fp);
    if (header == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: could not read BAM header in chromosome thread\n";
        sam_close(sam_fp);
        return chromResults;
    }

    // Load BAM index for random access
    // If explicit index URL is provided, use it; otherwise auto-detect
    hts_idx_t* idx = NULL;
    if (!opt::bamindexurl.empty()) {
        idx = sam_index_load2(sam_fp, bampath.c_str(), opt::bamindexurl.c_str());
    } else {
        idx = sam_index_load(sam_fp, bampath.c_str());
    }

    if (idx == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Warning: could not load BAM index for " << bampath << "\n";
        if (!opt::bamindexurl.empty()) {
            std::cerr << "Warning: tried explicit index URL: " << opt::bamindexurl << "\n";
        }
        sam_hdr_destroy(header);
        sam_close(sam_fp);
        return chromResults;
    }

    // Create iterator for this chromosome (or unmapped reads if tid == -1)
    hts_itr_t* iter = NULL;
    if (tid == -1) {
        // Query unmapped reads using HTS_IDX_NOCOOR
        iter = sam_itr_queryi(idx, HTS_IDX_NOCOOR, 0, 0);
    } else {
        // Query specific chromosome
        iter = sam_itr_queryi(idx, tid, 0, INT32_MAX);
    }

    if (iter == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: could not create iterator for chromosome " << chr_name << std::endl;
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(sam_fp);
        return chromResults;
    }

    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[chromosome " << chr_name << "] Starting processing\n";
    }

    // Allocate BAM alignment structure
    bam1_t* b = bam_init1();
    int nprocessed = 0;

    // Iterate through reads in this chromosome
    int ret;
    while ((ret = sam_itr_next(sam_fp, iter, b)) >= 0) {
        // Extract read group tag
        std::string tag = opt::unknown;
        if (rggroups) {
            uint8_t* rg_tag = bam_aux_get(b, "RG");
            if (rg_tag != NULL) {
                tag = std::string(bam_aux2Z(rg_tag));
            } else {
                continue; // skip reads without RG tag
            }
        }

        // skip reads with readgroup not defined in result template
        if (chromResults.find(tag) == chromResults.end()) {
            continue;
        }

        // for exome, exclude reads mapped to the exome regions.
        if (isExome) {
            range rg;
            rg.first = b->core.pos;
            rg.second = b->core.pos + b->core.l_qseq;
            std::string chrm = refID2Name(b->core.tid);

            if (chrm != "-1") {
                std::map<std::string, std::vector<range> >::iterator chrmit = opt::exomebed.find(chrm);
                if (chrmit == opt::exomebed.end()) {
                    chromResults[tag].n_exreadsChrUnmatched += 1;
                } else {
                    std::vector<range>::iterator itend = opt::exomebed[chrm].end();
                    std::map<std::string, std::vector<range>::iterator>::iterator lastfoundchrmit = lastfound.find(chrm);
                    if (lastfoundchrmit == lastfound.end()) {
                        lastfound[chrm] = chrmit->second.begin();
                    }
                    searchhint = lastfound[chrm];
                    std::vector<range>::iterator itsearch = searchRange(searchhint, itend, rg);
                    if (itsearch != itend) {
                        searchhint = itsearch;
                        chromResults[tag].n_exreadsExcluded += 1;
                        lastfound[chrm] = searchhint;
                        continue;
                    }
                }
            }
        }

        // Extract sequence
        std::string queryBases;
        uint8_t* seq = bam_get_seq(b);
        for (int i = 0; i < b->core.l_qseq; i++) {
            queryBases += seq_nt16_str[bam_seqi(seq, i)];
        }

        // Extract flags
        bool isMapped = !(b->core.flag & BAM_FUNMAP);
        bool isDuplicate = (b->core.flag & BAM_FDUP);

        // Update counts
        chromResults[tag].numTotal += 1;

        if (isMapped) {
            chromResults[tag].numMapped += 1;
        }

        if (isDuplicate) {
            chromResults[tag].numDuplicates += 1;
        }

        double gc = calcGC(queryBases);
        int ptn_count = countMotif(queryBases, opt::PATTERN, opt::PATTERN_REV);

        // Dynamically resize telcounts vector if needed for longer reads
        if (static_cast<size_t>(ptn_count) >= chromResults[tag].telcounts.size()) {
            chromResults[tag].telcounts.resize(ptn_count + 1, 0);
        }
        chromResults[tag].telcounts[ptn_count] += 1;

        if (gc >= ScanParameters::GC_LOWERBOUND && gc <= ScanParameters::GC_UPPERBOUND) {
            // get index for GC bin.
            int idx = floor((gc-ScanParameters::GC_LOWERBOUND)/ScanParameters::GC_BINSIZE);
            if (idx >= 0 && idx <= ScanParameters::GC_BIN_N-1 &&
                static_cast<size_t>(idx) < chromResults[tag].gccounts.size()) {
                chromResults[tag].gccounts[idx] += 1;
            }
        }

        nprocessed++;

        if (nprocessed % 5000000 == 0) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cerr << "[chromosome " << chr_name << "] processed " << nprocessed << " reads\n";
        }
    }

    // Clean up
    bam_destroy1(b);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(sam_fp);

    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[chromosome " << chr_name << "] Completed. Processed " << nprocessed << " reads\n";
    }

    return chromResults;
}

// Process a single BAM file with chromosome-level parallel processing (thread-safe function)
// Uses htslib for proper URL streaming support (S3, HTTP, HTTPS, FTP)
std::map<std::string, ScanResults> processSingleBam(const std::string& bampath, bool isExome) {

  // storing results for each read group (RG tag). use
  // read group ID as key.
  std::map<std::string, ScanResults> resultmap;
  // store where the overlap was last found in the case of exome seq
  std::map<std::string, std::vector<range>::iterator> lastfound;
  std::vector<range>::iterator searchhint;

  {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Start analysing BAM " << bampath << "\n";
  }

  // Open BAM file using htslib (supports URLs)
  samFile* sam_fp = sam_open(bampath.c_str(), "r");
  if (sam_fp == NULL) {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Error: could not open BAM file: " << bampath << std::endl;
    return resultmap;
  }

  // Read header
  sam_hdr_t* header = sam_hdr_read(sam_fp);
  if (header == NULL) {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Error: could not read BAM header from: " << bampath << std::endl;
    sam_close(sam_fp);
    return resultmap;
  }

  bool rggroups = false;

  // Parse read groups from header
  if(opt::ignorerg){ // ignore read groups
	  {
        std::lock_guard<std::mutex> lock(output_mutex);
	    std::cerr << "Treat all reads in BAM as if they were from a same sample" << std::endl;
      }
	  ScanResults results;
	  results.sample = opt::unknown;
	  resultmap[opt::unknown]=results;
  }else{
    std::map <std::string, std::string> readgroups;
    std::map <std::string, std::string> readlibs;

    // Parse @RG lines from header text
    const char* header_text = sam_hdr_str(header);
    if (header_text) {
      std::string hdr(header_text);
      free((void*)header_text);

      // Simple parser for @RG lines
      std::istringstream iss(hdr);
      std::string line;
      while (std::getline(iss, line)) {
        if (line.substr(0, 3) == "@RG") {
          rggroups = true;
          std::string rgid, sample, library;

          // Parse RG line fields
          std::istringstream linestream(line);
          std::string field;
          while (linestream >> field) {
            if (field.substr(0, 3) == "ID:") {
              rgid = field.substr(3);
            } else if (field.substr(0, 3) == "SM:") {
              sample = field.substr(3);
            } else if (field.substr(0, 3) == "LB:") {
              library = field.substr(3);
            }
          }

          if (!rgid.empty()) {
            readgroups[rgid] = sample.empty() ? opt::unknown : sample;
            readlibs[rgid] = library.empty() ? opt::unknown : library;
          }
        }
      }
    }

    if(rggroups){
	  {
        std::lock_guard<std::mutex> lock(output_mutex);
	    std::cerr<<"Specified BAM has "<< readgroups.size()<< " read groups" << std::endl;
      }

	  for(std::map<std::string, std::string>::iterator it = readgroups.begin(); it != readgroups.end(); ++it){
		  ScanResults results;
		  std::string rgid = it -> first;
		  results.sample = it -> second;
		  results.lib = readlibs[rgid];
		  resultmap[rgid]=results; //results are identified by RG tag.
	  }

	}else{
	  {
        std::lock_guard<std::mutex> lock(output_mutex);
	    std::cerr << "Warning: can't find RG tag in the BAM header" << std::endl;
	    std::cerr << "Warning: treat all reads in BAM as if they were from a same sample" << std::endl;
      }
	  ScanResults results;
	  results.sample = opt::unknown;
	  results.lib = opt::unknown;
	  resultmap[opt::unknown]=results;
	}
  }

  // Always use parallel chromosome processing (indexed BAM required)
  unsigned int threads_per_chrom = opt::threads_per_file;

  // Save n_targets and chromosome names
  int n_targets = header->n_targets;
  std::vector<std::string> chr_names;
  chr_names.reserve(n_targets);
  for (int i = 0; i < n_targets; i++) {
    chr_names.push_back(std::string(header->target_name[i]));
  }

  // Check if BAM index is available (required for parallel processing)
  hts_idx_t* test_idx = NULL;
  if (!opt::bamindexurl.empty()) {
    test_idx = sam_index_load2(sam_fp, bampath.c_str(), opt::bamindexurl.c_str());
  } else {
    test_idx = sam_index_load(sam_fp, bampath.c_str());
  }

  if (test_idx == NULL) {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Error: BAM index (.bai) not found for " << bampath << "\n";
    if (!opt::bamindexurl.empty()) {
      std::cerr << "Error: Explicit index URL was provided: " << opt::bamindexurl << "\n";
    }
    std::cerr << "Error: BAM index is required. Please index your BAM file using 'samtools index'.\n";
    std::cerr << "Error: For S3/HTTP URLs, use -i/--bam-index to specify the index URL.\n";
    sam_hdr_destroy(header);
    sam_close(sam_fp);
    return resultmap;
  }

  {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "BAM index found. Using " << threads_per_chrom << " threads for chromosome-level parallel processing\n";
  }

  // Parallel chromosome processing
  {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Found " << n_targets << " chromosomes/contigs in BAM header\n";
    std::cerr << "Processing chromosomes in parallel with " << threads_per_chrom << " threads\n";
  }

    // Process chromosomes in batches to limit concurrent threads
    std::vector<std::future<std::map<std::string, ScanResults>>> chromFutures;
    chromFutures.reserve(threads_per_chrom);

    // Process unmapped reads first (tid = -1)
    {
      std::lock_guard<std::mutex> lock(output_mutex);
      if (opt::unmappedonly) {
        std::cerr << "Processing ONLY unmapped reads (unmapped-only mode)...\n";
      } else {
        std::cerr << "Processing unmapped reads first...\n";
      }
    }
    chromFutures.push_back(std::async(std::launch::async,
        processChromosome, bampath, -1, "*", isExome, rggroups, resultmap));

    // Wait for unmapped reads to complete
    auto unmappedResult = chromFutures.front().get();
    for (const auto& pair : unmappedResult) {
      const std::string& rg = pair.first;
      if (resultmap.find(rg) != resultmap.end()) {
        add_results(resultmap[rg], pair.second);
      }
    }
    chromFutures.clear();

    // If unmapped-only mode, skip chromosome processing
    if (opt::unmappedonly) {
      // Don't clean up - child threads have their own copies and htslib may share internal structures
      std::lock_guard<std::mutex> lock(output_mutex);
      std::cerr << "Unmapped-only mode: Skipping mapped chromosome processing.\n";
      std::cerr << "Completed scanning BAM (unmapped reads only)\n";
      return resultmap;
    }

    // Now process mapped chromosomes in parallel
    int tid = 0;
    while (tid < n_targets) {
      // Launch up to threads_per_chrom tasks
      while (chromFutures.size() < threads_per_chrom && tid < n_targets) {
        std::string chr_name = chr_names[tid];  // Use saved chromosome name
        chromFutures.push_back(std::async(std::launch::async,
            processChromosome, bampath, tid, chr_name, isExome, rggroups, resultmap));
        tid++;
      }

      // Wait for at least one task to complete before launching more
      if (chromFutures.size() >= threads_per_chrom && tid < n_targets) {
        auto chromResult = chromFutures.front().get();
        for (const auto& pair : chromResult) {
          const std::string& rg = pair.first;
          if (resultmap.find(rg) != resultmap.end()) {
            add_results(resultmap[rg], pair.second);
          }
        }
        chromFutures.erase(chromFutures.begin());
      }
    }

    // Collect results from remaining chromosome threads ONE AT A TIME
    // to avoid memory issues with many large maps
    {
      std::lock_guard<std::mutex> lock(output_mutex);
      std::cerr << "[scan] Aggregating results from " << chromFutures.size()
                << " remaining chromosome threads...\n";
    }

    size_t future_idx = 0;
    for (auto& future : chromFutures) {
      {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[scan] Aggregating future " << (future_idx + 1) << " of " << chromFutures.size() << "...\n";
      }

      auto chromResult = future.get();

      // Aggregate results immediately
      for (const auto& pair : chromResult) {
        const std::string& rg = pair.first;
        if (resultmap.find(rg) != resultmap.end()) {
          add_results(resultmap[rg], pair.second);
        }
      }

      // chromResult will be destroyed here, freeing memory before next iteration
      future_idx++;
    }

    // Clear the futures vector to free any remaining resources
    chromFutures.clear();

  // Don't clean up main thread's header/sam_fp/test_idx - child threads have their own copies
  // and htslib may share internal structures. Let OS clean up on process exit.

  // Parallel processing complete
  {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "Completed scanning BAM\n";
  }
  return resultmap;
}

int scanBam()
{

  std::vector< std::map<std::string, ScanResults> > resultlist;
  bool isExome = opt::exomebedfile.size()==0? false: true;

  std::cout << opt::bamlist << "\n";
  std::cout << opt::bamlist.size() << " BAMs" <<  std::endl;

  // Check for URLs in the BAM list
  bool has_urls = false;
  for (const auto& bam : opt::bamlist) {
    if (isURL(bam)) {
      has_urls = true;
      std::cerr << "Detected URL: " << bam << "\n";
    }
  }

  // Determine the number of parallel BAM files to process
  unsigned int parallel_files = opt::parallel_files;
  unsigned int max_threads = std::thread::hardware_concurrency();

  // Auto-detect if parallel_files is 0
  if (parallel_files == 0) {
    parallel_files = max_threads > 0 ? max_threads : 1;
  }

  // Validate parallel_files
  if (parallel_files < 1) {
    parallel_files = 1;
  }
  if (max_threads > 0 && parallel_files > max_threads) {
    std::cerr << "Warning: Requested " << parallel_files << " parallel files, but only "
              << max_threads << " hardware threads available. Using " << max_threads << ".\n";
    parallel_files = max_threads;
  }

  std::cerr << "Processing up to " << parallel_files << " BAM files in parallel\n";
  if (opt::threads_per_file > 1) {
    std::cerr << "Using " << opt::threads_per_file << " threads per BAM file for read processing\n";
  }

  if (has_urls) {
    std::cerr << "Note: BAM files will be streamed from URLs. BamTools supports http://, https://, and ftp:// URLs.\n";
  }

  // Process BAM files in parallel using thread pool
  if (parallel_files == 1 || opt::bamlist.size() == 1) {
    // Single-threaded mode for file processing
    for(std::size_t i=0; i<opt::bamlist.size(); i++) {
      std::map<std::string, ScanResults> resultmap = processSingleBam(opt::bamlist[i], isExome);
      resultlist.push_back(resultmap);
    }
  } else {
    // Multi-threaded mode for file processing
    resultlist.resize(opt::bamlist.size());
    std::vector<std::future<std::map<std::string, ScanResults>>> futures;

    for(std::size_t i=0; i<opt::bamlist.size(); i++) {
      // Wait if we have too many files being processed
      while(futures.size() >= parallel_files) {
        for(auto it = futures.begin(); it != futures.end(); ) {
          if(it->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready) {
            it = futures.erase(it);
          } else {
            ++it;
          }
        }
      }

      // Launch new thread to process this BAM file
      futures.push_back(std::async(std::launch::async, processSingleBam, opt::bamlist[i], isExome));
    }

    // Collect results in order
    size_t result_idx = 0;
    for(auto& future : futures) {
      resultlist[result_idx++] = future.get();
    }
  }

  if(opt::onebam){
    merge_results_by_readgroup(resultlist);
  }

  outputresults(resultlist);

  if(isExome){
    printlog(resultlist);
  }

  std::cerr << "Completed writing results\n";

  return 0;
}

// Legacy single-threaded version removed - now using htslib with URL support

void printlog(std::vector< std::map<std::string, ScanResults> > resultlist){

	for(size_t i =0; i< resultlist.size(); i++){
		const auto& rmap = resultlist[i];
		for(std::map<std::string, ScanResults>::const_iterator it= rmap.begin();
				it != rmap.end(); ++it){

			const std::string& rg = it ->first;
			const ScanResults& result = it -> second;
			std::cout << "BAM:" << rg << std::endl;
			std::cout << "	chr ID unmatched reads: " << result.n_exreadsChrUnmatched << std::endl;
			std::cout << "	exome reads excluded: " << result.n_exreadsExcluded << std::endl;
		}
	}

}


void printout(const std::string& rg, const ScanResults& result, std::ostream* pWriter, size_t max_telcounts_size){

	*pWriter << rg << ScanParameters::FIELD_SEP;
	*pWriter << result.lib << ScanParameters::FIELD_SEP;
	*pWriter << result.sample << ScanParameters::FIELD_SEP;
	*pWriter << result.numTotal << ScanParameters::FIELD_SEP;
	*pWriter << result.numMapped << ScanParameters::FIELD_SEP;
	*pWriter << result.numDuplicates << ScanParameters::FIELD_SEP;

	double telLenEstimate = calcTelLength(result);
	if(telLenEstimate==-1){
		std::cerr << "Telomere length estimate unknown. No read was found with telomeric GC composition.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(telLenEstimate>1000000){
		std::cerr << "Telomere length estimate unknown. Possibly due to not enough representation of genome.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(telLenEstimate==0){
		std::cerr << "Telomere length estimate unknown. No read contains more than " << opt::tel_k << " telomere repeats.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}
	else{
		*pWriter << telLenEstimate << ScanParameters::FIELD_SEP;
	}

	// Output ALL telomere counts found, padding with 0 if this result has fewer than max
	for (size_t j = 0; j < max_telcounts_size; ++j){
		if (j < result.telcounts.size()) {
			*pWriter << result.telcounts[j] << ScanParameters::FIELD_SEP;
		} else {
			*pWriter << 0 << ScanParameters::FIELD_SEP;
		}
	}

	for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
		*pWriter << result.gccounts[k] << ScanParameters::FIELD_SEP;
	}
	*pWriter << "\n";
}


int outputresults(std::vector< std::map<std::string, ScanResults> > resultlist){

	// Find maximum telcounts size across all results to determine how many TEL columns to output
	size_t max_telcounts_size = ScanParameters::TEL_MOTIF_N; // Default minimum
	for (const auto& resultmap : resultlist) {
		for (const auto& pair : resultmap) {
			if (pair.second.telcounts.size() > max_telcounts_size) {
				max_telcounts_size = pair.second.telcounts.size();
			}
		}
	}

	std::ostream* pWriter;
	bool tostdout = opt::outputfile.empty() ? true:false;
	if(tostdout){
		pWriter = &std::cout;
	}else{
		pWriter = createWriter(opt::outputfile);
	}

	if(opt::writerheader){
		// Generate dynamic headers based on maximum telcounts found
		*pWriter << ScanParameters::LABEL_RG << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_LB << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_SAMPLE << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_TOTAL << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_MAPPED << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_DUP << ScanParameters::FIELD_SEP;
		*pWriter << ScanParameters::LABEL_LEN << ScanParameters::FIELD_SEP;

		// Output TEL headers dynamically based on actual max found
		for(size_t i=0; i < max_telcounts_size; i++){
			*pWriter << ScanParameters::LABEL_TEL << i << ScanParameters::FIELD_SEP;
		}

		// Output GC headers (fixed number)
		for(int i=0; i<ScanParameters::GC_BIN_N; i++){
			*pWriter << ScanParameters::LABEL_GC << i << ScanParameters::FIELD_SEP;
		}
		*pWriter << "\n";
	}

	ScanResults mergedrs;
	std::string grpnames = "";

	for(size_t i=0; i < resultlist.size();++i){

		const std::map<std::string, ScanResults>& resultmap = resultlist[i];

		// if merge read groups, take weighted average for all measures
		bool domg = opt::mergerg && resultmap.size() > 1? true:false;

		for(std::map<std::string, ScanResults>::const_iterator it= resultmap.begin();
				it != resultmap.end(); ++it){

			const std::string& rg = it ->first;
			const ScanResults& result = it -> second;

			if(domg){
				if(grpnames.size()==0){
					grpnames += rg;
				}else{
					grpnames += "|"+rg;
				}

				mergedrs.sample = result.sample;
				mergedrs.numTotal += result.numTotal;
				mergedrs.numMapped += result.numMapped * result.numTotal;
				mergedrs.numDuplicates += result.numDuplicates * result.numTotal;
				mergedrs.telLenEstimate += calcTelLength(result)* result.numTotal;

				// Ensure mergedrs vectors are large enough before accessing
				if (result.telcounts.size() > mergedrs.telcounts.size()) {
					mergedrs.telcounts.resize(result.telcounts.size(), 0);
				}
				if (result.gccounts.size() > mergedrs.gccounts.size()) {
					mergedrs.gccounts.resize(result.gccounts.size(), 0);
				}

				for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
					mergedrs.telcounts[j] += result.telcounts[j]* result.numTotal;
				}
				for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
					mergedrs.gccounts[k] += result.gccounts[k]* result.numTotal;
				}
				continue;
			}else{
				printout(rg, result, pWriter, max_telcounts_size);
			}
		}

		//in this case calculate weighted average
		if(domg){

			mergedrs.numMapped /= mergedrs.numTotal;
			mergedrs.numDuplicates /= mergedrs.numTotal;
			mergedrs.telLenEstimate /= mergedrs.numTotal;

			for (std::size_t j = 0, max = mergedrs.telcounts.size(); j != max; ++j){
				mergedrs.telcounts[j] /= mergedrs.numTotal;
			}
			for (std::size_t k = 0, max = mergedrs.gccounts.size(); k != max; ++k){
				mergedrs.gccounts[k] /= mergedrs.numTotal;
			}

			mergedrs.numTotal  /= resultmap.size();

			printout(grpnames, mergedrs, pWriter, max_telcounts_size);
		};
	}

	if(!tostdout){
		delete pWriter;
	}
	return 0;
}

double calcTelLength(const ScanResults& results){

	float acc = 0;
	for(std::size_t i=opt::tel_k, max = results.telcounts.size(); i !=max; ++i){
		acc += results.telcounts[i];
	}

	float gc_tel = 0;
	for(std::size_t i=0, max = results.gccounts.size(); i !=max; ++i){
		float gc1=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*i;
		float gc2=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*(i+1);
		if(gc1 >= ScanParameters::GC_TELOMERIC_LOWERBOUND && gc2 <= ScanParameters::GC_TELOMERIC_UPPERBOUND ){
			gc_tel += results.gccounts[i];
		}
	}

	if(gc_tel == 0){
		return -1;
	}

	return (acc/gc_tel)*float(ScanParameters::GENOME_LENGTH_AT_TEL_GC)/ScanParameters::LENGTH_UNIT/ScanParameters::TELOMERE_ENDS;
}



int countMotif(const std::string& read, const std::string& pattern, const std::string& pattern_revcomp){

	int motifcount1 = 0;
	int motifcount2 = 0;

	// Cache pattern sizes to avoid repeated calls
	const size_t pattern_len = pattern.size();
	const size_t revcomp_len = pattern_revcomp.size();

	size_t p1 = read.find(pattern, 0);
	while(p1 != std::string::npos)
	{
	    p1 = read.find(pattern, p1 + pattern_len);
	    motifcount1++;
	}

	size_t p2 = read.find(pattern_revcomp, 0);
	while(p2 != std::string::npos)
	{
	    p2 = read.find(pattern_revcomp, p2 + revcomp_len);
	    motifcount2++;
	}

	return motifcount1 > motifcount2 ? motifcount1 : motifcount2;

}

double calcGC(const std::string& seq)
{
    const size_t len = seq.size();
    if (len == 0) return 0.0;  // Handle empty sequence

    double num_gc = 0.0;
    for(char c : seq)
    {
        if(c == 'C' || c == 'G')
            ++num_gc;
    }
    return num_gc / len;
}

void update_pattern(){

	opt::PATTERN = ScanParameters::PATTERN;
	opt::PATTERN_REV = ScanParameters::PATTERN_REVCOMP;

    // update total motif counts when read length and/or pattern have been specified by user
    ScanParameters::TEL_MOTIF_N = ScanParameters::READ_LENGTH/ScanParameters::PATTERN.size() +1;
	if(opt::tel_k < 1 or opt::tel_k > ScanParameters::TEL_MOTIF_N-1){
		std::cerr << "k out of bound. k must be an integer from 1 to " <<  ScanParameters::TEL_MOTIF_N-1 << "\n";
		exit(EXIT_FAILURE);
	}

}


//
// Handle command line arguments
//
void parseScanOptions(int argc, char** argv)
{

	std::string bamlistfile =  "";
	std::string rev = "";

	Headers hd;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'f':
            	arg >> bamlistfile; break;
            case 'i':
            	arg >> opt::bamindexurl; break;
            case 'o':
            	arg >> opt::outputfile; break;
            case 'H':
            	opt::writerheader=false; break;
            case 'm':
                opt::mergerg = true; break;
            case 'u':
                opt::ignorerg = true; break;
            case 'U':
                opt::unmappedonly = true; break;
            case 'w':
                opt::onebam = true; break;
            case 'h':
        		for(size_t h=0; h<hd.headers.size();h++){
        			std::cout << hd.headers[h] << ScanParameters::FIELD_SEP;
        		}
        		std::cout << "\n";
        		exit(EXIT_SUCCESS);
            case 'k':
            	arg >> opt::tel_k;
            	break;
            case 'j':
            	arg >> opt::parallel_files;
            	break;
            case 'T':
            	arg >> opt::threads_per_file;
            	if(opt::threads_per_file < 1){
            		std::cerr << "Threads per file must be at least 1\n";
            		exit(EXIT_FAILURE);
            	}
            	break;
            case 'r':
				arg >> ScanParameters::READ_LENGTH;
				if(ScanParameters::READ_LENGTH <= 0 || ScanParameters::READ_LENGTH > 100000){
					std::cerr << "please specify valid read length that is greater than 0 and length than 100kb" << "\n";
					exit(EXIT_FAILURE);
				}
				break;
            case 'p':

				break;
            case 'z':
            	arg >> ScanParameters::PATTERN;
				ScanParameters::PATTERN_REVCOMP = reverseComplement(ScanParameters::PATTERN);
				std::cerr << "use user specified pattern " <<  ScanParameters::PATTERN << "\n";
				std::cerr << "reverse complement " <<  ScanParameters::PATTERN_REVCOMP << "\n";
            	break;
            case 'e':
            	arg >> opt::exomebedfile;
            	opt::exomebed = readBedAsVector(opt::exomebedfile);
            	std::cout << "loaded "<< opt::exomebed.size() << " exome regions \n"<< std::endl;
//            	std::cout << opt::exomebed << "\n";
            	break;

            case OPT_HELP:
                std::cout << TELSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TELSEQ_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    update_pattern();

    // deal with cases of API usage:
    // telseq a.bam b.bam c.bam ...
    // | telseq
    // telseq -f

    if (argc - optind < 1) // no argument specified
    {
    	// check if it is from pipe
		if(!isatty(fileno(stdin))){
			std::string line;
			while (std::getline(std::cin, line))
			{
				if(line.empty()){
					continue;
				}
				opt::bamlist.push_back(line);
			}
		}else if(bamlistfile.empty() ){ // check if not from a pipe, -f must be spceified
//    		std::cerr << SUBPROGRAM ": No BAM specified. Please specify BAM either directly, by using -f option or piping BAM file path.\n";
    		std::cout << TELSEQ_USAGE_MESSAGE;
    		exit(EXIT_FAILURE);
    	}
    }

    else if (argc - optind >= 1) // if arguments are specified
    {
    	// -f has higher priority, when specified, ignore arguments.
    	if(bamlistfile.empty()){
    		for(int i = optind; i < argc; ++i ){
    		    opt::bamlist.push_back(argv[i]);
    		}
    	}
    }

    // read in bamlist
    if(!bamlistfile.empty()){
        size_t filesize = getFilesize(bamlistfile);
        if(filesize == 0)
        {
            std::cerr << PROGRAM_BIN ": BAMLIST file specified by -f is empty\n";
            exit(EXIT_FAILURE);
        }

        std::istream* pReader = createReader(bamlistfile);
        std::string line;

        while(getline(*pReader, line))
        {
        	if(line.empty()){
        		continue;
        	}
            opt::bamlist.push_back(line);
        }

        size_t bamsize = opt::bamlist.size();
        if(bamsize == 0 ){
            std::cerr << PROGRAM_BIN ": Could not find any sample in BAMLIST file specified.\n";
            exit(EXIT_FAILURE);
        }
        delete pReader;
    }
}



int main(int argc, char** argv)
{
	Timer* pTimer = new Timer("scan BAM");
	parseScanOptions(argc, argv);
    scanBam();
    delete pTimer;
    return 0;
}
