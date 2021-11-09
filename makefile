# Varialble Definitions and Flags:
CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Compiler is mpicc    NOTE: you need to have an MPI implementation installed on your computer
LD = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Linker is mpicc also

CFALGS = -Wall -Wunused-value -fcommon
CDEBUG = -g -backtrace

MPI_COMPILE_FLAGS = $(bash mpicc --showme:compile)
MPI_LINK_FLAGS = $(bash mpicc --showme:link)

#LNETCFD = -lnetcdf

# Linker flags go here. Currently there aren't any, but if we'll switch to
# code optimization, we might add "-s" here to strip debug info and symbols.
LDFLAGS =
BOUND_FLAG =

# List of source codes and generated object files:
SRCS = \
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
	./src/MESH.c \
	./src/NRUTIL.c \
	./src/PARABOLA.c \
	./src/PRINT.c \
	./src/READ_INPUT.c \
	./src/READ_TOPOGRAPHY.c \
	./src/SURFACES.c \
	./src/TOPOfromTXT.c \
	./src/TOPO_USER_FUNCTION.c \
	./src/WRITE_OUTPUT.c

#	./src/DOMAIN_DECOMP.c \
	./src/MINMAXVAL.c \
	./src/elliptic_solver.c \
	./src/mympi_init.c \

OBJS = \
	./src/main.o \
	./src/ALMOST_EQUAL.o \
	./src/BUILD_CONN.o \
	./src/BUILD_GRID.o \
	./src/BUILD_LGL.o \
	./src/GAUSSJ.o \
	./src/GRID2CONN.o \
	./src/GRID_COORD.o \
	./src/INTERPOLATE.o \
	./src/LINSPACE.o \
	./src/MEMORY.o \
	./src/MESH.o \
	./src/NRUTIL.o \
	./src/PARABOLA.o \
	./src/PRINT.o \
	./src/READ_INPUT.o \
	./src/READ_TOPOGRAPHY.o \
	./src/SURFACES.o \
	./src/TOPOfromTXT.o \
	./src/TOPO_USER_FUNCTION.o \
	./src/WRITE_OUTPUT.o


BIN=./bin
SRC=./src
OBJ=./obj

# ProgRam executable file name:
EXE = $(BIN)/MMesh3D-V2.a

# Top-level rule, to compile everything
all: $(EXE)

# Rule to Link the program
$(EXE): $(OBJS)
	$(LD) $(CFLAGS) $(MPI_LINK_FLAGS) -L/usr/include $(OBJS) -o $(EXE)

# Rule for the file "main.o":
main.o: main.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/main.c

clean:
	rm -rf  $(BIN)/* *~ $(SRC)/*.o
	clear

