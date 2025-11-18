#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TelSeq], [0.1.0], [igor@sequencing.com])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Telseq/telseq.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for header files.
AC_LANG([C++])

AC_ARG_WITH(htslib,
        AS_HELP_STRING([--with-htslib=DIR], [specify directory containing the htslib library]),
        [ if test -d "$withval" ; then
            htslib_ldflags="-Wl,-rpath,$withval/lib -L$withval/lib"
    		htslib_include="-I$withval/include"
		else
          	AC_MSG_ERROR([Cannot stat directory $withval])
      	fi ]
        )

# Set compiler flags.
#AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror -Wno-unknown-pragmas")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Wno-unknown-pragmas -g -O0")
#AC_SUBST(CXXFLAGS, "-O3")
#AC_SUBST(CFLAGS, "-O3")

AC_SUBST(CXXFLAGS, "-std=c++11")
AC_SUBST(CFLAGS, "-std=c++11")

AC_SUBST(CPPFLAGS, "$CPPFLAGS $htslib_include $gzstream_include")
AC_SUBST(LDFLAGS, "$htslib_ldflags $gzstream_ldflags $LDFLAGS")

#AC_MSG_NOTICE(["CPPFLAGS is $CPPFLAGS"])
#AC_MSG_NOTICE(["LDFLAGS is $LDFLAGS"])

# We always need to specify to link in htslib
AC_SUBST(LIBS, "$LIBS -lhts -lz")

#AC_MSG_NOTICE(["LIBS is $LIBS"])

# Make sure the htslib headers can be found
AC_CHECK_HEADER([htslib/sam.h],,[AC_MSG_ERROR(
[The htslib library must be installed. \
You can specify its path with the --with-htslib=PATH option])])

AC_CONFIG_FILES([Makefile
		Util/Makefile
    Telseq/Makefile])

AC_OUTPUT
