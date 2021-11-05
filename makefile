# Varialble Definitions and Flags:
CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Compiler is mpicc    NOTE: you need to have an MPI implementation installed on your compute
CXX = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpic++ #Compiler is mpicc    NOTE: you need to have an MPI implementation installed on your computer
#LD = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Linker is mpicc als
LD = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpic++ #Linker is mpicc also

CFALGS = -Wall -Wunused-value -fcommon
CDEBUG = -g -backtrace

MPI_COMPILE_FLAGS = $(bash mpic++ --showme:compile)
MPI_LINK_FLAGS = $(bash mpic++ --showme:link)

#LNETCFD = -lnetcdf

# Linker flags go here. Currently there aren't any, but if we'll switch to
# code optimization, we might add "-s" here to strip debug info and symbols.
LDFLAGS =
BOUND_FLAG =

# List of source codes and generated object files:
SRCS = \
	./src/main.cpp \
	./src/MYVECTOR.cpp \
	./src/MEMORY.cpp \
	./src/nrutil.c
#	./src/READ_INPUT.c	\
#	./src/PRINT.c		\
#	./src/BUILD_LGL.c


OBJS = \
	./src/main.o \
	./src/MYVECTOR.o \
	./src/MEMORY.o \
	./src/nrutil.o
#	./src/READ_INPUT.o	\
#	./src/PRINT.o		\
#	./src/BUILD_LGL.o

BIN=./bin
SRC=./src

# ProgRam executable file name:
EXE = $(BIN)/MMesh3D-V2++.a

# Top-level rule, to compile everything
all: $(EXE)

# Rule to Link the program
$(EXE): $(OBJS)
	$(CXX) $(CFLAGS) $(MPI_LINK_FLAGS) -L/usr/include $(OBJS) -o $(EXE)

clean:
	rm -rf $(SRC)/*.o $(SRC)/*~ $(BIN)/* *~
	clear

