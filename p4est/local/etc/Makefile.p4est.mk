
# This file is part of p4est
# Use `include /path/to/Makefile.p4est.mk' in your Makefile
# to use p4est in your project without autotools

prefix = /Users/simone/Work/Codes/mmesh3d/github/MMesh3D/p4est/local
exec_prefix = ${prefix}
p4est_sysconfdir = ${prefix}/etc

include ${p4est_sysconfdir}/Makefile.sc.mk

# P4EST_CC and P4EST_CFLAGS may not be necessary for your project
P4EST_CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc
P4EST_CFLAGS = -g -O2

# These pull in p4est but none of its dependencies
P4EST_PKG_CPPFLAGS = -I${prefix}/include
P4EST_PKG_LDFLAGS = -L${exec_prefix}/lib
P4EST_PKG_LIBS = -lp4est

# These pull in everything needed by p4est
P4EST_CPPFLAGS =  $(SC_PKG_CPPFLAGS) $(P4EST_PKG_CPPFLAGS)
P4EST_LDFLAGS =  $(SC_PKG_LDFLAGS) $(P4EST_PKG_LDFLAGS)
P4EST_LIBS = $(P4EST_PKG_LIBS) $(SC_PKG_LIBS) -lz 
