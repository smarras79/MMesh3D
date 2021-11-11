/* src/pre_config.h.  Generated from pre_config.h.in by configure.  */
/* src/pre_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* C compiler */
#define CC "/Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc"

/* C compiler flags */
#define CFLAGS "-g -O2"

/* C preprocessor */
#define CPP "/Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc -E"

/* C preprocessor flags */
#define CPPFLAGS ""

/* DEPRECATED (use SC_ENABLE_DEBUG instead) */
/* #undef DEBUG */

/* enable debug mode (assertions and extra checks) */
/* #undef ENABLE_DEBUG */

/* Undefine if: While the default alignment is sizeof (void *), this switch
   will choose the standard system malloc. For custom alignment use
   --enable-memalign=<bytes> */
#define ENABLE_MEMALIGN 1

/* Define to 1 if we are using MPI */
#define ENABLE_MPI 1

/* Define to 1 if we can use MPI_COMM_TYPE_SHARED */
#define ENABLE_MPICOMMSHARED 1

/* Define to 1 if we are using MPI I/O */
#define ENABLE_MPIIO 1

/* Define to 1 if we are using MPI_Init_thread */
#define ENABLE_MPITHREAD 1

/* Define to 1 if we can use MPI_Win_allocate_shared */
#define ENABLE_MPIWINSHARED 1

/* enable OpenMP: Using --enable-openmp without arguments does not specify any
   CFLAGS; to supply CFLAGS use --enable-openmp=<OPENMP_CFLAGS>. We check
   first for linking without any libraries and then with -lgomp; to avoid the
   latter, specify LIBS=<OPENMP_LIBS> on configure line */
/* #undef ENABLE_OPENMP */

/* enable POSIX threads: Using --enable-pthread without arguments does not
   specify any CFLAGS; to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
   We check first for linking without any libraries and then with -lpthread;
   to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line */
/* #undef ENABLE_PTHREAD */

/* Undefine if: disable non-thread-safe internal debug counters */
#define ENABLE_USE_COUNTERS 1

/* Undefine if: replace array/dmatrix resize with malloc/copy/free */
#define ENABLE_USE_REALLOC 1

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the `backtrace' function. */
#define HAVE_BACKTRACE 1

/* Define to 1 if you have the `backtrace_symbols' function. */
#define HAVE_BACKTRACE_SYMBOLS 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the `fsync' function. */
#define HAVE_FSYNC 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Have we found function pthread_create. */
/* #undef HAVE_LPTHREAD */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Have we found function omp_get_thread_num. */
/* #undef HAVE_OPENMP */

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the `qsort_r' function. */
#define HAVE_QSORT_R 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the `strtoll' function. */
#define HAVE_STRTOLL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Have we found function adler32_combine. */
#define HAVE_ZLIB 1

/* Define to 1 on a bigendian machine */
/* #undef IS_BIGENDIAN */

/* Linker flags */
#define LDFLAGS ""

/* Libraries */
#define LIBS "-lz "

/* minimal log priority */
/* #undef LOG_PRIORITY */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* DEPRECATED (use SC_ENABLE_MEMALIGN instead) */
#define MEMALIGN 1

/* desired alignment of allocations in bytes */
#define MEMALIGN_BYTES (SC_SIZEOF_VOID_P)

/* DEPRECATED (use SC_ENABLE_MPI instead) */
#define MPI 1

/* DEPRECATED (use SC_ENABLE_MPIIO instead) */
#define MPIIO 1

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* DEPRECATED (use SC_ENABLE_OPENMP instead) */
/* #undef OPENMP */

/* Name of package */
#define PACKAGE "libsc"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "info@p4est.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libsc"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libsc 2.8.1.44-954f"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libsc"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.8.1.44-954f"

/* DEPRECATED (use SC_WITH_PAPI instead) */
/* #undef PAPI */

/* Use builtin getopt */
/* #undef PROVIDE_GETOPT */

/* DEPRECATED (use SC_ENABLE_PTHREAD instead) */
/* #undef PTHREAD */

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* DEPRECATED (use SC_ENABLE_USE_COUNTERS instead) */
#define USE_COUNTERS 1

/* DEPRECATED (use SC_ENABLE_USE_REALLOC instead) */
#define USE_REALLOC 1

/* Version number of package */
#define VERSION "2.8.1.44-954f"

/* Package major version */
#define VERSION_MAJOR 2

/* Package minor version */
#define VERSION_MINOR 8

/* Package point version */
#define VERSION_POINT 1.44-954f

/* enable Flop counting with papi */
/* #undef WITH_PAPI */

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */
