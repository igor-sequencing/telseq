# TelSeq - Enhanced Fork with S3/HTTP Streaming Support

TelSeq is a software that estimates telomere length from whole genome sequencing data (BAMs).

This is an enhanced fork with the following improvements:
- **S3/HTTP URL streaming support** with explicit BAM index URLs
- **Chromosome-level parallelization** for faster processing
- **Fixed critical parallel processing bug** that caused incorrect read counts
- Uses **htslib** instead of BamTools for better URL support and compatibility

## Original Project

This fork is based on the original TelSeq by Zhihao Ding:
- Original repository: [https://github.com/zd1/telseq](https://github.com/zd1/telseq)
- This fork: [https://github.com/igor-sequencing/telseq](https://github.com/igor-sequencing/telseq)

### Citation

_Estimating telomere length from whole genome sequence data_

Zhihao Ding; Massimo Mangino; Abraham Aviv; Tim Spector; Richard Durbin.
Nucleic Acids Research 2014; doi: 10.1093/nar/gku181
[http://nar.oxfordjournals.org/content/42/9/e75](http://nar.oxfordjournals.org/content/42/9/e75)

## What's New in This Fork

### Bug Fixes
- **Fixed parallel processing bug**: In the original implementation, when using multi-threading, read counts were incorrectly multiplied by the number of threads due to accumulated results being passed to each thread.

### New Features
- **S3/HTTP URL streaming**: Process BAM files directly from S3 or HTTP URLs without downloading
- **Explicit BAM index URLs**: Use `-i/--bam-index` to specify the .bai file URL separately (required for presigned URLs)
- **Chromosome-level parallelization**: Use `-T/--threads-per-file` to process chromosomes in parallel within each BAM file
- **File-level parallelization**: Use `-j/--jobs` to process multiple BAM files in parallel
- **Unmapped-only mode**: Use `-U` to quickly process only unmapped reads and exit

### Migration Notes
- The deprecated `-t/--threads` option has been removed
- Use `-j/--jobs` for file-level parallelism (processing multiple BAMs concurrently)
- Use `-T/--threads-per-file` for read-level parallelism (processing chromosomes within a BAM concurrently)

## Dependencies

### Required
- **htslib** (version 1.10 or later) - replaces the BamTools dependency
  ```bash
  # macOS
  brew install htslib

  # Ubuntu/Debian
  apt-get install libhts-dev

  # From source
  git clone https://github.com/samtools/htslib.git
  cd htslib
  make && make install
  ```

- **Modern C++ compiler** (GCC 4.8+ or Clang with C++11 support)
  ```bash
  # macOS
  brew install gcc

  # Ubuntu/Debian
  apt-get install build-essential
  ```

- **Autotools** (automake, autoconf)
  ```bash
  # macOS
  brew install automake autoconf

  # Ubuntu/Debian
  apt-get install automake autoconf
  ```

## Compilation

### Quick Start

```bash
cd src
./autogen.sh
./configure
make
```

The executable binary will be at `src/Telseq/telseq`.

### Custom htslib Location

If htslib is installed in a non-standard location:

```bash
./configure --with-htslib=/path/to/htslib
make
```

The `/path/to/htslib` directory should contain `lib` and `include` subdirectories.

### Specify Compiler

If you need to use a specific compiler:

```bash
export CXX=/path/to/g++
export CC=/path/to/gcc
./configure
make
```

## Usage

### Basic Usage

```bash
# Analyze local BAM files
telseq sample1.bam sample2.bam

# Analyze BAMs listed in a file
telseq -f bamlist.txt

# From pipe
cat bamlist.txt | telseq
```

### S3/HTTP URL Streaming (NEW)

Process BAM files directly from URLs without downloading:

```bash
# With auto-detected index (for public URLs with standard .bai naming)
telseq https://example.com/sample.bam

# With explicit index URL (required for presigned URLs)
telseq -i https://example.com/sample.bam.bai https://example.com/sample.bam

# S3 presigned URLs (BAM and BAI must be specified separately)
telseq -i "https://bucket.s3.region.amazonaws.com/sample.bam.bai?..." \
       "https://bucket.s3.region.amazonaws.com/sample.bam?..."
```

### Parallel Processing (NEW)

```bash
# Process multiple BAM files in parallel (file-level parallelism)
telseq -j 4 sample1.bam sample2.bam sample3.bam sample4.bam

# Process chromosomes in parallel within each BAM (read-level parallelism)
# Requires BAM index (.bai file)
telseq -T 8 sample.bam

# Combine both: process 2 BAMs concurrently, using 4 threads per BAM
telseq -j 2 -T 4 sample1.bam sample2.bam

# Auto-detect number of cores
telseq -j 0 sample1.bam sample2.bam  # -j 0 uses all available cores
```

### Quick Unmapped Read Analysis (NEW)

Process only unmapped reads for quick telomere length estimation:

```bash
# Requires BAM index
telseq -U sample.bam

# Works with URLs too
telseq -U -i "https://example.com/sample.bam.bai" "https://example.com/sample.bam"
```

### Output Options

```bash
# Write to file
telseq -o output.tsv sample.bam

# Or redirect stdout
telseq sample.bam > output.tsv 2> log.txt

# Suppress header (useful when combining multiple outputs)
telseq -H sample1.bam sample2.bam > combined.tsv

# Print header only
telseq -h > header.tsv
```

### Read Group Handling

```bash
# Ignore read groups (treat all reads as one sample)
telseq -u sample.bam

# Merge read groups (weighted average)
telseq -m sample.bam
```

## Command Line Options

### Input/Output
- `-f, --bamlist=FILE` - File containing list of BAM paths (one per line)
- `-o, --output-dir=FILE` - Output file (default: stdout)
- `-i, --bam-index=URL` - URL to BAM index file (.bai), required for S3/HTTP URLs with -T > 1

### Parallelization (NEW)
- `-j, --jobs=INT` - Max number of BAM files to process in parallel (0 = auto-detect cores, default: 0)
- `-T, --threads-per-file=INT` - Number of threads for processing reads within each BAM file (requires BAM index, default: 1)

### Analysis Options
- `-U` - Process ONLY unmapped reads and exit (requires BAM index, useful for quick analysis)
- `-u` - Ignore read groups (treat all reads as one sample)
- `-m` - Merge read groups (weighted average)
- `-k INT` - Threshold for telomeric repeats (default: 7)
- `-e, --exomebed=FILE` - Exclude exome regions in BED format

### Output Formatting
- `-H` - Remove header line
- `-h` - Print header line only
- `--help` - Display full help message
- `--version` - Display version information

### Advanced Options
- `-r INT` - Read length (default: 100)
- `-z PATTERN` - Use custom pattern for searching (e.g., TTAGGG)
- `-w` - Consider all BAMs as one single BAM

## Output Format

| Column | Definition |
|--------|------------|
| ReadGroup | Read group ID from BAM header (RG tag) |
| Library | Sequencing library from BAM header (LB tag) |
| Sample | Sample name from BAM header (SM tag) |
| Total | Total number of reads |
| Mapped | Number of mapped reads (SAM flag 0x4) |
| Duplicates | Number of duplicate reads (SAM flag 0x400) |
| LENGTH_ESTIMATE | Estimated telomere length (in kb) |
| TEL0-TEL16 | Read counts by number of TTAGGG/CCCTAA repeats |
| GC0-GC9 | Read counts by GC content (40%-60% in 2% bins) |

## Examples

### Example 1: Local BAM Files

```bash
# Simple analysis
telseq sample.bam

# Multiple files with parallel processing
telseq -j 4 -T 2 sample1.bam sample2.bam sample3.bam sample4.bam
```

### Example 2: S3 Streaming

```bash
# Public S3 bucket
telseq https://s3.amazonaws.com/bucket/sample.bam

# Presigned S3 URLs (common in production environments)
BAM_URL="https://bucket.s3.region.amazonaws.com/sample.bam?AWSAccessKeyId=...&Signature=..."
BAI_URL="https://bucket.s3.region.amazonaws.com/sample.bam.bai?AWSAccessKeyId=...&Signature=..."
telseq -i "$BAI_URL" "$BAM_URL" -T 4
```

### Example 3: Quick Analysis with Unmapped Reads

```bash
# Fast telomere length estimation using only unmapped reads
telseq -U -u sample.bam

# With S3 streaming
telseq -U -u -i "$BAI_URL" "$BAM_URL"
```

### Example 4: Batch Processing

```bash
# Create BAM list
ls /data/*.bam > bamlist.txt

# Process all BAMs in parallel
telseq -j 0 -f bamlist.txt -o results.tsv
```

## Performance Tips

1. **Use parallel processing**: `-j` and `-T` can significantly speed up processing
   - `-j` for multiple files (I/O bound)
   - `-T` for large files (CPU bound)

2. **For quick estimates**: Use `-U` to process only unmapped reads (~10-100x faster)

3. **S3 streaming**:
   - Always use `-i` to specify the index URL explicitly for presigned URLs
   - Use `-T` for parallel chromosome processing when streaming from S3
   - Consider network bandwidth when choosing thread count

4. **Memory usage**: Higher thread counts use more memory (each thread opens the BAM independently)

## Docker

### Build Docker Image

```bash
docker build -t telseq-docker github.com/igor-sequencing/telseq
```

### Run with Docker

```bash
# Local BAM file
docker run -v /path/to/bam:/data telseq-docker /data/sample.bam

# With URL (no volume mount needed)
docker run telseq-docker https://example.com/sample.bam
```

## Troubleshooting

### BAM Index Not Found
```
Warning: BAM index (.bai) not found
```
**Solution**: When using URLs, explicitly specify the index with `-i`:
```bash
telseq -i https://example.com/sample.bam.bai https://example.com/sample.bam
```

### Incorrect Read Counts (Old Bug - Fixed)
If using an older version of telseq with multi-threading, read counts may be multiplied by the number of threads. **This bug has been fixed in this fork.**

### Compilation Errors
```
fatal error: htslib/sam.h: No such file or directory
```
**Solution**: Install htslib or specify its location:
```bash
./configure --with-htslib=/path/to/htslib
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues at:
[https://github.com/igor-sequencing/telseq](https://github.com/igor-sequencing/telseq)

## License

See LICENSE file for details.

## Contact

For issues specific to this fork:
- GitHub Issues: [https://github.com/igor-sequencing/telseq/issues](https://github.com/igor-sequencing/telseq/issues)

For the original TelSeq project:
- Original author: zhihao.ding at gmail.com
- Original repository: [https://github.com/zd1/telseq](https://github.com/zd1/telseq)
