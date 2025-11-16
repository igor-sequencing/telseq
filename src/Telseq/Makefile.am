bin_PROGRAMS = telseq

telseq_CPPFLAGS = \
	-I$(top_srcdir)/Util

telseq_CXXFLAGS = \
	-pthread

telseq_LDFLAGS = \
	-pthread

telseq_LDADD = \
	$(top_builddir)/Util/libutil.a \
	-lhts

telseq_SOURCES = telseq.cpp telseq.h 