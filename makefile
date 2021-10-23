# Varialble Definitions and Flags:
CC = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Compiler is mpicc    NOTE: you need to have an MPI implementation installed on your computer
LD = /Users/simone/mylibs/openmpi-4.1.0/build_gcc9/bin/mpicc #Linker is mpicc also

CFALGS = -Wall -Wunused-value -fcommon
CDEBUG = -g

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
	./src/MEMORY.c \
	./src/READ_INPUT.c \
	./src/READ_TOPOGRAPHY.c \
	./src/BUILD_GRID.c \
	./src/BUILD_CONN.c \
	./src/GRID_COORD.c \
	./src/PRINT.c \
	./src/GRID2CONN.c \
	./src/SURFACES.c \
	./src/nrutil.c \
	./src/linspace.c \
	./src/parabola.c \
	./src/gaussj.c \
	./src/TOPOfromTXT.c \
	./src/INTERPOLATE.c \
	./src/minmaxval.c \
	./src/topo_user_function.c \
	./src/BUILD_LGL.c \
	./src/WRITE_OUTPUT.c \
	./src/elliptic_solver.c \
	./src/PSUP.c \
	./src/mympi_init.c \
	./src/DOMAIN_DECOMP.c

OBJS = \
	./src/main.o \
	./src/MEMORY.o \
	./src/READ_INPUT.o \
	./src/READ_TOPOGRAPHY.o \
	./src/BUILD_GRID.o \
	./src/BUILD_CONN.o \
	./src/GRID_COORD.o \
	./src/PRINT.o \
	./src/GRID2CONN.o \
	./src/SURFACES.o \
	./src/nrutil.o \
	./src/linspace.o \
	./src/parabola.o \
	./src/gaussj.o \
	./src/TOPOfromTXT.o \
	./src/INTERPOLATE.o \
	./src/minmaxval.o \
	./src/topo_user_function.o \
	./src/BUILD_LGL.o \
	./src/WRITE_OUTPUT.o \
	./src/elliptic_solver.o \
	./src/PSUP.o \
	./src/mympi_init.o \
	./src/DOMAIN_DECOMP.o

BIN=./bin

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

# Rule for the file "READ_INPUT.o":
READ_INPUT.o: READ_INPUT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/READ_INPUT.c

# Rule for the file "READ_TOPOGRAPHY.o":
READ_TOPOGRAPHY.o: READ_TOPOGRAPHY.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/READ_TOPOGRAPHY.c

# Rule for the file "READ_INPUT.o":
MEMORY.o: MEMORY.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/MEMORY.c

# Rule for the file "BUILD_GRID.o":
BUILD_GRID.o: BUILD_GRID.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/BUILD_GRID.c

# Rule for the file "BUILD_CONN.o":
BUILD_CONN.o: BUILD_CONN.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/BUILD_CONN.c

# Rule for the file "GRID_COORD.o":
GRID_COORD.o: GRID_COORD.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/GRID_COORD.c

# Rule for the file "PRINT.o":
PRINT.o: PRINT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/PRINT.c

# Rule for the file "GRID2CONN.o":
GRID2CONN.o: GRID2CONN.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/GRID2CONN.c

# Rule for the file "SURFACES.o":
SURFACES.o: SURFACES.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/SURFACES.c

# Rule for the file "WRITE_OUTPUT.o":
WRITE_OUTPUT.o: WRITE_OUTPUT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/WRITE_OUTPUT.c

# Rule for the file "nrutil.o"
nrutil.o: nrutil.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/nrutil.c

# Rule for the file "linspace.o"
linspace.o: linspace.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/linspace.c

# Rule for the	 file "parabola.o"
parabola.o: parabola.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/parabola.c

# Rule for the	 file "gaussj.o"
gaussj.o: gaussj.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/gaussj.c

# Rule for the	 file "TOPOfromTXT.o"
TOPOfromTXT.o: TOPOfromTXT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/TOPOfromTXT.c

# Rule for the	 file "INTERPOLATE.o"
INTERPOLATE.o: INTERPOLATE.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/INTERPOLATE.c

# Rule for the	 file "minmaxval.o"
minmaxval.o: minmaxval.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/minmaxval.c

# Rule for the   file "topo_user_function.o"
topo_user_function.o: topo_user_function.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/topo_user_function.c

# Rule for the file "visit_writer.o"
BUILD_LGL.o: BUILD_LGL.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/BUILD_LGL.c

# Rule for the file "visit_writer.o"
elliptic_solver.o: elliptic_solver.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/elliptic_solver.c

# Rule for the file "PSUP.o"                                                                                 
PSUP.o: PSUP.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/PSUP.c

# Rule for the file "visit_writer.o"
visit_writer.o: visit_writer.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/visit_writer.c

# Rule for the file "mympi_init.o"
mympi_init.o: mympi_init.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/mympi_init.c

# Rule for the file "DOMAIN_DECOMP.o"
DOMAIN_DECOMP.o: DOMAIN_DECOMP.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c ./src/DOMAIN_DECOMP.c

destroy:
	rm -rf $(OBJ) $(BIN)/*

clean:
	rm -rf $(OBJ) $(BIN)/* *~ *.o *.msh *.dom.dat
	clear

