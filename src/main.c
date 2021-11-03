#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

int main(int argc, char** argv) {
    // Initialize the MPI environment
    //MPI_Init(NULL, NULL);
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
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

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

	/*************************************************************************************
	 *READ INPUT FILE
	 *************************************************************************************/
	READ_INPUT(inputfile);

	PRINT_INFO();

	/*************************************************************************************
	 * READ THE TOPOGRAPHY FILE or BUILD THE BOTTOM BOUNDARY BY A USER-DEFINED FUNCTION:
	 *************************************************************************************/
	READ_TOPOGRAPHY();
	double L;
	int ix;
	int nx;
	double *x;

	st_legendre Legendre;
	
	FILE *file_id;
	nx = 100;
	file_id = fopen("LEGENDRE.dat", "w");       
	x = dvector(1,nx);
	dlinspace(-1, 1, nx, x, "n");
	
	for (ix=1; ix<=nx; ix++) {

	    Legendre = LegendreAndDerivative(nop, x[ix]);
	    fprintf(file_id, "% - .8f %.8f %f\n", x[ix], Legendre.legendre, Legendre.dlegendre);
	}
	free_dvector(x, 1, nx);
  	fclose(file_id);
	
	st_lgl lgl;
	lgl.size    = ngl;
	lgl.coords  = (double*) malloc( sizeof(double) * lgl.size);
	lgl.weights = (double*) malloc( sizeof(double) * lgl.size);

	LegendreGaussLobattoNodesAndWeights(lgl, nop);
	//LegendreGaussNodesAndWeights(lgl, nop);
	
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
	
    }

     if (irank == 0) {
	 /*****************************************************
	  * Free memory
	  *****************************************************/
	 MEMORY_DEALLOCATE(1);
	 MEMORY_DEALLOCATE(2);
	 MEMORY_DEALLOCATE(0);
     }
    
    
    // Finalize the MPI environment.
    MPI_Finalize();
}
