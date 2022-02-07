/*
#ifndef P4_TO_P8
#define P4_TO_P8
#endif
*/
#include <stdbool.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GLOBAL_VARS.h" //GLOBAL VARIABLES
#include "MYDEFINE.h"    //GLOBAL CONSTANTS
#include "MYSTRUCTS.h"
#include "MESH.h"

//Functions headers
#include "ALMOST_EQUAL.h"
#include "BUILD_CONN.h"
#include "BUILD_GRID.h"
#include "BUILD_LGL.h"
#include "GMSH_IO.h"
#include "GRID_COORD.h"
#include "GRID2CONN.h"
#include "INTERPOLATE.h"
#include "MEMORY.h"
#include "NRUTIL.h"
#include "PRINT.h"
#include "READ_INPUT.h"
#include "READ_TOPOGRAPHY.h"
#include "SURFACES.h"
#include "WRITE_OUTPUT.h"
#include "TOPOfromTXT.h"
#include "TOPO_USER_FUNCTION.h"

//#include "BUILD_GRID_SPHERE.h"

//#include "linspace.h"
//#include "parabola.h"
//#include "GAUSSJ.h"
//#include "MINMAXVAL.h"
//#include "elliptic_solver.h"
//#include "mympi_init.h"
//#include "DOMAIN_DECOMP.h"
//#include "visit_writer.h"


int main(int argc, char** argv) {

    st_lgl lgl;
    /*************************************************************************************
     *Initialize MPI and get the process' rank
     *************************************************************************************/
    MPI_Init(&argc, &argv);
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    // Get the name of the processor
    //char processor_name[MPI_MAX_PROCESSOR_NAME];
    //int name_len;
    //MPI_Get_processor_name(processor_name, &name_len);

    if (irank == 0) {

	/*****************************************************************************
	 * Command-line arguments:
	 *****************************************************************************/
   	if(argc == 1){
	    //Remember that the first argument of argc is the name of the program by default, 
	    //thus argc == 1 correswponds to that and the program recognizes that you did not put any input
	    printf("\n You did not enter the input file\n");
	    printf(" Call the program as:\n");
	    printf("\n ./meshgen.a [inputfile.inp][alya]\n\n");
	    printf(" where\n");
	    printf("           [inputfile.inp] is the file of input parameters\n");
	    printf("           [alya] is an OPTIONAL argument indicating whether\n");
	    printf("                  you want to or not to generate the alya-readable file\n\n");
	    return 1;
	}
	else{

	    PRINT_WELCOME_MESSAGE();
	    
	    /*****************************************************
	     * Allocate string pointers and vectors used in input:
	     *****************************************************/
	    MEMORY_ALLOCATE(0);
	    
	    if(argc == 2 ){
		strcpy(inputfile, argv[1]);
	    }
	    else if(argc == 3){
		strcpy(inputfile, argv[1]);
		strcpy(print_alya, argv[2]);
	    }
	}
	
	/*************************************************************************************
	 *READ INPUT FILE
	 *************************************************************************************/
	READ_INPUT(inputfile);
	PRINT_INFO();

	/*************************************************************************************
	 *COMPUTE HIGH-ORDER NODES and WEIGHTS. Stored in struct: lgl.coords and lgl.weights
	 *************************************************************************************/	
	if (lread_external_grid == 0){
	    
	    /*************************************************************************************
	     * READ THE TOPOGRAPHY FILE or BUILD THE BOTTOM BOUNDARY BY A USER-DEFINED FUNCTION:
	     *************************************************************************************/
	    READ_TOPOGRAPHY();
	    
	    /*************************************************************************************
	     * Dynamic memory allocation of COORDS and CONN related arrays
	     *************************************************************************************/
	    MEMORY_ALLOCATE(1);
	    MEMORY_ALLOCATE(10);
	    
	    /*************************************************************************************
	     * Build GRID coordinates and connectivity matrix:
	     *************************************************************************************/
	    BUILD_GRID();
	    
	    /*************************************************************************************
	     * Build the Connectivity matrix
	     *************************************************************************************/
	    BUILD_CONN();

	    //apply_smoothing();
	    
	    /*************************************************************************************
	     * Split the initial domain into 'mpiprocs' processors:
	     *************************************************************************************/
	    ////DOMAIN_DECOMP(rank);
	    
	    //Add high order nodes:
	    if (nop > 1) {
		printf(" !!! ERROR in main:!!!\n");
		printf(" !!! NOP > 1 for native mesh generator is not working yet.\n");
		printf(" !!! Still missing CONN_BDY_FACES from the grid constraction.\n");
		printf(" !!! The program will EXIT now.\n");
		printf(" !!! Set `nop 1` in the input file `%s`\n\n", inputfile);
		exit(-1);
		/*
		//CGNS_ORDERING(CONN, nelem);
		BUILD_EDGES(CONN, nelem);
		
		MEMORY_ALLOCATE(11);
		ADD_HIGH_ORDER_NODES();
	    	
		MEMORY_DEALLOCATE(6);
		*/
	    }
	    
	} else {
	    GMSH_IO(external_grid_file_name);

	    //apply_smoothing();
	    
	    //CGNS_ORDERING(CONN, nelem);
	    BUILD_EDGES(CONN, nelem);
	    
	    MEMORY_ALLOCATE(11);
	    ADD_HIGH_ORDER_NODES();
	    	    
	    MEMORY_DEALLOCATE(6);

	} //END reading external grid
	
	/*************************************************************************************
	 * Write output to file (VTK, ALYA, etc.)
	 *************************************************************************************/
	WRITE_OUTPUT(irank);	   
	
	/*****************************************************
	 * Free memory
	 *****************************************************/
	printf(" #------------------------------------------------------------------#\n");
	MEMORY_DEALLOCATE(1);
	MEMORY_DEALLOCATE(10);
	MEMORY_DEALLOCATE(2);
	MEMORY_DEALLOCATE(0);
	if (lread_external_grid == 1){
	    MEMORY_DEALLOCATE(11);    
	    MEMORY_DEALLOCATE(7);
	    MEMORY_DEALLOCATE(5);
	}
	printf(" #------------------------------------------------------------------#\n");
    }
    
    // Finalize the MPI environment.
    MPI_Finalize();
}

