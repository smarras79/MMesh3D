/* 
 * This functions allocate/deallocate 
 * the memory for the global pointers/arrays
 * when needed.
 *
 * Simone Marras, Aug 2012
 */
#include <errno.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NRUTIL.h"
#include "MYDEFINE.h"
#include "GLOBAL_VARS.h"

/*
 * Allocate memory
 */
int MEMORY_ALLOCATE(int flag)
{

  int inode, iel;

  if(flag == 0)
    {
      inputfile  = (char*) malloc(32 * sizeof(char *));
      print_alya = (char*) malloc(3 * sizeof(char *));
      print_vtk  = (char*) malloc(3 * sizeof(char *));
      print_gmsh  = (char*) malloc(3 * sizeof(char *));
      fname      = (char*) malloc(23 * sizeof(char));
      input_inp  = (char*) malloc(8 * sizeof(char));
  
      outfile_msh_vtk  = (char*) malloc(96 * sizeof(char *));
      outfile_msh_gmsh  = (char*) malloc(96 * sizeof(char *));
      outfile_msh_alya = (char*) malloc(96 * sizeof(char *));
  
      parameters = dvector(0,NUMBER_OF_PARAMETERS);

      //Allocation of the strings of the above lines:
      for(i=0; i<PROB_ENTRIES; i++)
	problem[i] = (char *) malloc(32 * sizeof(char *));
      
    }
  
  if(flag == 1)
    {
      printf(" EL_NODES %d\n", EL_NODES);
      /*coordinates*/
      COORDS   = dmatrix(0,nnodes,0,nsd+1);
      COORDS1d = dvector(0,nnodes*(nsd+1));
      
      /*coonectivity*/
      CONN = imatrix(0,nelem,0,EL_NODES+1); //For hexa only for now
      for(iel=0; iel<=nelem; iel++){
	for(inode=0; inode<=EL_NODES+1; inode++){
	  CONN[iel][inode] = 0;
	}
      }
      
      /*Element type*/
      ELTYPE  = ivector(0,nelem+1);
  
      /*Bdy flag*/
      BDYFLAG = dmatrix(0, nnodes,0,nsd+1);
      
    }

  if(flag == 2)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      BDY_COORDS = dmatrix(0,NBDY_NODES,0,nsd);
    }

  if(flag == 3)
    {
      rarray = dvector(0,nnodes);
    }

  if(flag == 4)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      x = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);
      y = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);
      z = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);
    }
  
  return 0;
}


/*
 * Free memory
 */
int MEMORY_DEALLOCATE(int flag)
{
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
      
      free_dvector(parameters,0,NUMBER_OF_PARAMETERS);
      free_dvector(COORDS1d,  0,nnodes*(nsd+1));
      free_ivector(ELTYPE,    0,nelem+1);
      free_dmatrix(COORDS,    0,nnodes,     0,nsd+1);
      free_imatrix(CONN,      0,nelem,      0,EL_NODES+1);
      free_dmatrix(BDY_COORDS,0,NBDY_NODES, 0,nsd);
    }
  
  if(flag == 3)
    free_dvector(rarray, 0,nnodesx*nnodesy);
  
  if(flag == 4)
    {     
      //Allocate space for the array of boundary coordinates to be read in input:
      free_d3tensor(x,0,nnodesx,0,nnodesy,0,nnodesz);
      free_d3tensor(y,0,nnodesx,0,nnodesy,0,nnodesz);
      free_d3tensor(z,0,nnodesx,0,nnodesy,0,nnodesz);
    }

return 0;    
}

