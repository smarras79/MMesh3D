/* 
 * This functions allocate/deallocate 
 * the memory for the global pointers/arrays
 * when needed.
 *
 * Simone Marras, Aug 2012
 */

//C++
#include <iostream>
#include "MEMORY.hpp"
#include "MYVECTOR.hpp"

//C
#ifdef __cplusplus
extern "C" {
#endif
#include "global_vars.h"
#ifdef __cplusplus
}
#endif


/*
 * Allocate memory
 */
int MEMORY_ALLOCATE(int flag, ...)
{
  int inode, iel;
  
  /* if(flag == 0)
    {
      inputfile  = new char[32];
      print_alya = new char[3];
      print_vtk  = new char[3];
      print_gmsh = new char[3];
      fname      = new char[23];
      input_inp  = new char[8];
  
      outfile_msh_vtk  = new char[96];
      outfile_msh_gmsh = new char[96];
      outfile_msh_alya = new char[96];

      parameters = new double[NUMBER_OF_PARAMETERS + 1];
      
      //Allocation of the strings of the above lines:
      for(i=0; i<PROB_ENTRIES; i++)
	  problem[i] = new char[PROB_ENTRIES];
      
	  }*/
  
  if(flag == 1)
    {
      printf(" sEL_NODES %d\n", EL_NODES);
      nnodes = 11;
      
      //coordinates
      size_t nrows = nnodes + 1;
      size_t ncols = nsd + 1;
      
      myVector lCOORDS(nrows*ncols);
      myMatrix MCOORDS(nrows, ncols);

      lCOORDS.Print();
      
      MCOORDS.Print();
    }   
      /*
      COORDS    = new double*[nrows]; //allocate double* indexcing array
      COORDS[0] = new double[nrows*ncols];
      for(int i=1; i<=nrows; i++)
	  COORDS[i] = COORDS[0] + i*ncols;

      //connectivity
      nrows = nelem + 1;
      ncols = EL_NODES + 1;

      CONN    = new double*[nrows]; //allocate double* indexcing array
      CONN[0] = new double[nrows*ncols];
      for(int i=1; i<=nrows; i++)
	  CONN[i] = CONN[0] + i*ncols;
      
      //COORDS   = dmatrix(1,nnodes, 1,nsd+1);
      COORDS1d = new double[nnodes*(nsd + 1) + 1];
            
      //Element type
      ELTYPE  = new int[nelem+1];
  
      //Bdy flag
      nrows = nnodes + 1;
      ncols = nsd + 1;

      BDYFLAG    = new double*[nrows]; //allocate double* indexcing array
      BDYFLAG[0] = new double[nrows*ncols];
      for(int i=1; i<=nrows; i++)
	  BDYFLAG[i] = BDYFLAG[0] + i*ncols;
      
    }
  /*
  if(flag == 2)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      BDY_COORDS = dmatrix(1,NBDY_NODES, 1,nsd);
    }

  if(flag == 3)
    {
      rarray = dvector(0,nnodes);
    }

  if(flag == 4)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      x = d3tensor(1,nnodesx, 1,nnodesy, 1,nnodesz);
      y = d3tensor(1,nnodesx, 1,nnodesy, 1,nnodesz);
      z = d3tensor(1,nnodesx, 1,nnodesy, 1,nnodesz);
    }
  */ 
  return 0;
}




/*
 * Free memory
 */
int MEMORY_DEALLOCATE(int flag, ...)
{
    /*
  if(flag == 0)
    {
      free(inputfile);
      free(print_alya);
      free(print_vtk);
      free(fname);
      free(input_inp);
  
      free(outfile_msh_vtk);
      free(outfile_msh_alya);
  
      for(i=0; i<PROB_ENTRIES; i++)
	free(problem[i]);
	}
  if(flag == 1)
    {
      
      free_dvector(parameters,1,NUMBER_OF_PARAMETERS);
      free_dmatrix(COORDS,   1,nnodes, 1,nsd+1);
      free_dvector(COORDS1d, 1,nnodes*(nsd+1));
      free_imatrix(CONN,   1,nelem,1,EL_NODES+1);
      free_ivector(ELTYPE, 0,nelem+1);
      free_dmatrix(BDY_COORDS, 1,NBDY_NODES, 1,nsd);
    }
  
  if(flag == 3)
    free_dvector(rarray, 0,nnodesx*nnodesy);
  
  if(flag == 4)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      free_d3tensor(x, 1,nnodesx, 1,nnodesy, 1,nnodesz);
      free_d3tensor(y, 1,nnodesx, 1,nnodesy, 1,nnodesz);
      free_d3tensor(z, 1,nnodesx, 1,nnodesy, 1,nnodesz);
    }
  */
return 0;    
}
