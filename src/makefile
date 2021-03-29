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
	main.c \
	MEMORY.c \
	READ_INPUT.c \
	READ_TOPOGRAPHY.c \
	BUILD_GRID.c \
	BUILD_CONN.c \
	GRID_COORD.c \
	PRINT.c \
	GRID2CONN.c \
	SURFACES.c \
	nrutil.c \
	linspace.c \
	parabola.c \
	gaussj.c \
	TOPOfromTXT.c \
	INTERPOLATE.c \
	minmaxval.c \
	topo_user_function.c \
	BUILD_LGL.c \
	WRITE_OUTPUT.c \
	elliptic_solver.c \
	PSUP.c \
	mympi_init.c \
	DOMAIN_DECOMP.c

OBJS = \
	main.o \
	MEMORY.o \
	READ_INPUT.o \
	READ_TOPOGRAPHY.o \
	BUILD_GRID.o \
	BUILD_CONN.o \
	GRID_COORD.o \
	PRINT.o \
	GRID2CONN.o \
	SURFACES.o \
	nrutil.o \
	linspace.o \
	parabola.o \
	gaussj.o \
	TOPOfromTXT.o \
	INTERPOLATE.o \
	minmaxval.o \
	topo_user_function.o \
	BUILD_LGL.o \
	WRITE_OUTPUT.o \
	elliptic_solver.o \
	PSUP.o \
	mympi_init.o \
	DOMAIN_DECOMP.o

# ProgRam executable file name:
EXE = MMesh3D-V2.a

# Top-level rule, to compile everything
all: $(EXE)

# Rule to Link the program
$(EXE): $(OBJS)
	$(LD) $(CFLAGS) $(MPI_LINK_FLAGS) -L/usr/include $(OBJS) -o $(EXE)

# Rule for the file "main.o":
main.o: main.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c main.c

# Rule for the file "READ_INPUT.o":
READ_INPUT.o: READ_INPUT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c READ_INPUT.c

# Rule for the file "READ_TOPOGRAPHY.o":
READ_TOPOGRAPHY.o: READ_TOPOGRAPHY.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c READ_TOPOGRAPHY.c

# Rule for the file "READ_INPUT.o":
MEMORY.o: MEMORY.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c MEMORY.c

# Rule for the file "BUILD_GRID.o":
BUILD_GRID.o: BUILD_GRID.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c BUILD_GRID.c

# Rule for the file "BUILD_CONN.o":
BUILD_CONN.o: BUILD_CONN.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c BUILD_CONN.c

# Rule for the file "GRID_COORD.o":
GRID_COORD.o: GRID_COORD.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c GRID_COORD.c

# Rule for the file "PRINT.o":
PRINT.o: PRINT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c PRINT.c

# Rule for the file "GRID2CONN.o":
GRID2CONN.o: GRID2CONN.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c GRID2CONN.c

# Rule for the file "SURFACES.o":
SURFACES.o: SURFACES.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c SURFACES.c

# Rule for the file "WRITE_OUTPUT.o":
WRITE_OUTPUT.o: WRITE_OUTPUT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c WRITE_OUTPUT.c

# Rule for the file "nrutil.o"
nrutil.o: nrutil.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c nrutil.c

# Rule for the file "linspace.o"
linspace.o: linspace.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c linspace.c

# Rule for the	 file "parabola.o"
parabola.o: parabola.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c parabola.c

# Rule for the	 file "gaussj.o"
gaussj.o: gaussj.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c gaussj.c

# Rule for the	 file "TOPOfromTXT.o"
TOPOfromTXT.o: TOPOfromTXT.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c TOPOfromTXT.c

# Rule for the	 file "INTERPOLATE.o"
INTERPOLATE.o: INTERPOLATE.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c INTERPOLATE.c

# Rule for the	 file "minmaxval.o"
minmaxval.o: minmaxval.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c minmaxval.c

# Rule for the   file "topo_user_function.o"
topo_user_function.o: topo_user_function.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c topo_user_function.c

# Rule for the file "visit_writer.o"
BUILD_LGL.o: BUILD_LGL.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c BUILD_LGL.c

# Rule for the file "visit_writer.o"
elliptic_solver.o: elliptic_solver.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c elliptic_solver.c

# Rule for the file "PSUP.o"                                                                                 
PSUP.o: PSUP.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c PSUP.c

# Rule for the file "visit_writer.o"
visit_writer.o: visit_writer.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c visit_writer.c

# Rule for the file "mympi_init.o"
mympi_init.o: mympi_init.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c mympi_init.c

# Rule for the file "DOMAIN_DECOMP.o"
DOMAIN_DECOMP.o: DOMAIN_DECOMP.c
	$(CC) $(CFLAGS) $(MPI_COMPILE_FLAGS) $(CDEBUG) -I/usr/include -c DOMAIN_DECOMP.c

destroy:
	rm -f $(OBJ) $(BIN)

clean:
	rm -f $(OBJ) *~ *.o *.msh *.dom.dat
	clear

