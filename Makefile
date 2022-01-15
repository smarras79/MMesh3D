#
# Set USER_NAME and set your own CC if not defined yet
#
#USER_NAME = "sm_eddy"
USER_NAME = "sm_imac"
#USER_NAME = "sm_macair"

#
# ADD YOUR CC value below with your choice of USER_NAME
#
$(info USER_NAME=$(USER_NAME))
ifeq ($(USER_NAME),"sm_macair")
	CC = /opt/homebrew/Cellar/openmpi/4.1.2/bin/mpicc
else ifeq ($(USER_NAME),"sm_imac")
	CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc
else ifeq ($(USER_NAME),"sm_eddy")
	CC = /opt/intel/oneapi/mpi/2021.3.0//bin/mpicc
else
missing_cc_error:
	@echo " ERROR in Makefile!"
	@echo "    CC rule not defined for $(USER_NAME)!"
	@echo "    Edit Makefile and add set your CC."
	exit
endif

MPI_COMPILE_FLAGS = $(bash mpicc --showme:compile)
MPI_LINK_FLAGS = $(bash mpicc --showme:link)

# define any compile-time flags
CFLAGS = -Wall -g -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comment -Wno-unused-function
#CFLAGS = -Wall -g -fbacktrace -Wunused-value -fcommon

# define any directories containing header files other than /usr/include
#

P4EST_DIR := 
#P4EST_DIR := /Users/simone/Work/Codes/mmesh3d/github/MMesh3D/p4est/local

INCLUDES = $(MPI_COMPILE_FLAGS) -I/usr/include
#INCLUDES = $(MPI_COMPILE_FLAGS) -I/usr/include -isystem$(P4EST_DIR)/include 

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#LFLAGS = -L/usr/local/lib -L$(P4EST_DIR)/lib
LFLAGS = 

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lm
#LIBS = -lp4est -lsc -lz -lm


SRC =./src
BIN =./bin
EXE = $(BIN)/MMesh3D.a

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

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
.PHONY: depend clean

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(SRC)/*.o $(SRC)/*~ $(EXE)
	clear
