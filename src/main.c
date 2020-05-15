/*****************************************************************************
 *  This is a structured mesh generator for domain with topography.
 *  On the original version of mesh generation in a box,
 *  the option of spherical grids was added.
 *  The part of spherical grid was partially taken from F. X. Giraldo's
 *  atmospheric code NUMA (http://faculty.nps.edu/fxgirald/Homepage/Link_to_NUMA.html)
 *
 *
 *  It is the option of generating elliptically smoothed grids with 
 *  quasi orthogonality.
 *
 * V2.0: 
 * - more input options
 * - cleaner main
 *
 * Simone Marras, August 2012
 * simone.marras@gmail.com
 *
 * Simone Marras, April 2014
 *****************************************************************************/
#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

int main(int argc, char *argv[])
{
  /* Declare local variables */
  int irank;
  int rank;
  int iter;
  int nprocs;
  int tag;
  int dest;

  double a, atot, amean, buff, err;
  
  /* Declare MPI variables */
  MPI_Status status;
  MPI_Request send_request; 
  MPI_Request reicv_request; 
  
  /* Initialize and define
   *   some values of certain scalars 
   */
  tag  = 0;
  dest = 0;
  rank = 0;
  
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
  //READ INPUT FILE:
  READ_INPUT(inputfile);
  
  /*************************************************************************************
   *Initialize MPI and get the process' rank
   *************************************************************************************/
  MPI_Init(&argc, &argv);
  
  /* Get process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Get number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &mpiprocs);
  
  //compute a on each proc:
  a = 10.0; // + (double)rank;
  
  if(rank != 0)
    {
      MPI_Isend(&a, 1, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD, &send_request);

      //MPI_Send(&a, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
      //printf(" Rank %d sends a = %lf to 0\n", rank, a);
    }
  else
    {
      atot = a;
      
      for(irank = 1; irank < mpiprocs; irank++)
	{
	  MPI_Irecv(&a, 1, MPI_DOUBLE, irank, tag, MPI_COMM_WORLD, &reicv_request);
	  //MPI_Recv(&a, 1, MPI_DOUBLE, irank, tag, MPI_COMM_WORLD, &status);

	  atot = atot + a;
	  //printf(" Rank 0 receives a = %lf from rank %d. Atot is now %lf\n", a, irank, atot);
	}
    }
  
  /*************************************************************************************
   * SUMMARY OF INPUT ENTRIES:
   *************************************************************************************/
  if(rank == 0) 
    PRINT_INFO();

  //TEMP SPHERICAL GRID:
  //BUILD_GRID_SPHERE(); //DO NOT USE YET!!! STILL TO BE CODED. NOT WORKING YET!!!
  
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
  //DOMAIN_DECOMP(rank);
  
  /*************************************************************************************
   * Build the Connectivity matrix
   *************************************************************************************/
  BUILD_CONN();
  
  /*************************************************************************************
   * Write output to file (VTK, ALYA, etc.)
   *************************************************************************************/
  //apply_smoothing();
  WRITE_OUTPUT(rank);

  /*************************************************************************************
   * Finalize MPI
   *************************************************************************************/
  MPI_Finalize();

  /*****************************************************************************************
   * Print message of ending to screen:
   *****************************************************************************************/
  if(rank == 0){
    printf("\n # !!! ENJOY YOUR NEW MESH !!! \n");
  }

  
  /*************************************************************************************
   * FREE dynamic memory
   *************************************************************************************/
  MEMORY_DEALLOCATE(0);
  
  return 0; 
}
