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
	/*****************************************************
	 * Allocate string pointers and vectors used in input:
	 *****************************************************/
	MEMORY_ALLOCATE(0);
	
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
	    if(argc == 2 ){
		strcpy(inputfile, argv[1]);
	    }
	    else if(argc == 3){
		strcpy(inputfile, argv[1]);
		strcpy(print_alya, argv[2]);
	    }
	}
    
	/*************************************************************************************
	 *Some more allocation after input arguments (nsd):
	 *************************************************************************************/
	MEMORY_ALLOCATE(2);

	if (irank == 0) {
	    /*************************************************************************************
	     *READ INPUT FILE
	     *************************************************************************************/
	    READ_INPUT(inputfile);
	    PRINT_INFO();

	    /*************************************************************************************
	     *COMPUTE HIGH-ORDER NODES and WEIGHTS. Stored in struct: lgl.coords and lgl.weights
	     *************************************************************************************/
	    //lgl = BUILD_LGL(nop);
	    //BarycentricWeights(lgl, nop);

	    int lread_gmsh = 0;
	    if (lread_gmsh == 0) {

		/*************************************************************************************
		 * READ THE TOPOGRAPHY FILE or BUILD THE BOTTOM BOUNDARY BY A USER-DEFINED FUNCTION:
		 *************************************************************************************/
		READ_TOPOGRAPHY();

		/*************************************************************************************
		 * Dynamic memory allocation of U,V,P,G,F and coordinates on the grid
		 *************************************************************************************/
		MEMORY_ALLOCATE(1);
	
		/*************************************************************************************
		 * Build GRID coordinates and connectivity matrix:
		 *************************************************************************************/	
		BUILD_GRID();
  
		/*************************************************************************************
		 * Split the initial domain into 'mpiprocs' processors:
		 *************************************************************************************/
		////DOMAIN_DECOMP(rank);
  
		/*************************************************************************************
		 * Build the Connectivity matrix
		 *************************************************************************************/
		BUILD_CONN();

		 
		/*************************************************************************************
		 * Write output to file (VTK, ALYA, etc.)
		 *************************************************************************************/
		//apply_smoothing();
		WRITE_OUTPUT(irank);	   

		/*****************************************************
		 * Free memory
		 *****************************************************/
		MEMORY_DEALLOCATE(1);
	    
	    } else {
		/*
		 * READ external *.gmsh file
		 */
		//READ_GMSH();
	   
	    }

	    //SOLVER de calor VA AQUI
	    
	}
    }


    
    if (irank == 0) {
	
	/*****************************************************
	 * Free memory
	 *****************************************************/
	MEMORY_DEALLOCATE(2);
	MEMORY_DEALLOCATE(0);
    }
    
    
    // Finalize the MPI environment.
    MPI_Finalize();
}
