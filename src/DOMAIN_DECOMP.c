/*
 * This function contains all the necessary functions
 * that perform domain decomposition.
 *
 * Since MMesh3D is logically cartesian, the Cartesian 
 * domain decomposer of MPI will be used.
 * 
 * Simone Marras, April 2014
 *
 */
#include "myinclude.h"
#include "global_vars.h"
#include "mydefine.h"

int DOMAIN_DECOMP(int irank)
{
  /* Declare local variables */
  MPI_Comm GRID_COMM;
  
  int      dim_sizes[nsd];   /* # of processes along each dimension             */
  int      periodic[nsd]; /* Periodicity on (logical, 1) or off (logical, 0) */
  int      reorder;          /* Reorder the processes order                     */

  int      npx, npy, npz;
  int      my_grid_rank;

  int      coordinates[nsd];

  /* Initialize values */
  reorder = 0;
  periodic[0] = 0;
  periodic[1] = 0;
  periodic[2] = 0;

  npx = (int)pow(mpiprocs, 1.0/3.0);
  npy = (int)pow(mpiprocs, 1.0/3.0);
  npz = (int)pow(mpiprocs, 1.0/3.0);
  
  if(irank==0){
    printf(" Nprocs = %d\n", mpiprocs);
    printf(" nprocx = %d\n", npx);
    printf(" nprocy = %d\n", npy);
    printf(" nprocz = %d\n", npz);
  }
  
  dim_sizes[0] = npx;  dim_sizes[1] = npy;  dim_sizes[2] = npz;
  if(!strncmp(problem[3], "On", 2))
    periodic[0] = 1;
  if(!strncmp(problem[4], "On", 2))
    periodic[1] = 1;
  if(!strncmp(problem[8], "On", 2))
    periodic[2] = 1;
  
  /* MPI_CART */
  if(lMETIS == 0 || lCART == 1)
    {
      /* Start with the calls to MPI_CART: */
      MPI_Cart_create(MPI_COMM_WORLD, nsd, dim_sizes, periodic, reorder, &GRID_COMM);
      
      /* Get the coords of the current process */
      MPI_Cart_coords(GRID_COMM, irank, nsd, coordinates);
      printf(" Irank %d has coordinates: %d, %d, %d\n", irank, coordinates[0], coordinates[1], coordinates[2]);
      
      //MPI_Cart_rank(GRID_COMM, coordinates, &irank);
      
    }
  else
    {
      /* Other decomposition methods can be added here 
       * e.g. METIS, ...
       */
      if(lMETIS == 1)
	{
	  printf(" METIS not implemented yet. MPI_CART will be used by default\n" );
	  
	}
      
    }
  
  return 0;
}
