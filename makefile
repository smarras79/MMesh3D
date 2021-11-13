#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc

MPI_COMPILE_FLAGS = $(bash mpicc --showme:compile)
MPI_LINK_FLAGS = $(bash mpicc --showme:link)

# define any compile-time flags
CFLAGS = -Wall -g -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comment -Wno-unused-function
#CFLAGS = -Wall -g -fbacktrace -Wunused-value -fcommon

# define any directories containing header files other than /usr/include
#

#P4EST_DIR := 
P4EST_DIR := /Users/simone/Work/Codes/mmesh3d/github/MMesh3D/p4est/local

#INCLUDES = $(MPI_COMPILE_FLAGS) -I/usr/include
INCLUDES = $(MPI_COMPILE_FLAGS) -I/usr/include -isystem$(P4EST_DIR)/include 

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -L/usr/local/lib -L$(P4EST_DIR)/lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
#LIBS =
LIBS = -lp4est -lsc -lz

# define the C source files
SRCS =  \
	./src/main.c \
	./src/ALMOST_EQUAL.c \
	./src/BUILD_CONN.c \
	./src/BUILD_GRID.c \
	./src/BUILD_LGL.c \
	./src/GAUSSJ.c \
	./src/GRID2CONN.c \
	./src/GRID_COORD.c \
	./src/INTERPOLATE.c \
	./src/LINSPACE.c \
	./src/MEMORY.c \
	./src/NRUTIL.c \
	./src/PARABOLA.c \
	./src/PRINT.c \
	./src/READ_INPUT.c \
	./src/READ_TOPOGRAPHY.c \
	./src/SURFACES.c \
	./src/TOPOfromTXT.c \
	./src/TOPO_USER_FUNCTION.c \
	./src/WRITE_OUTPUT.c \
	./src/GMSH_IO.c

#	./src/P4EST_API.c
#	./src/MESH.c \

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file
SRC=./src
BIN=./bin
EXE = $(BIN)/MMesh3D.a

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
.PHONY: depend clean

all: $(EXE)
	@echo "---------------------------------------------"
	@echo  COMPILATION of $(BIN)/$(EXE) was SUCCESSFUL. 
	@echo "---------------------------------------------"

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(SRC)/*.o *~ $(EXE)
	clear

depend: $(SRCS)
	makedepend $(INCLUDES) $^

