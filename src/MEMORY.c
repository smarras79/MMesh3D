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

#include "GLOBAL_VARS.h"
#include "MESH.h"
#include "MYDEFINE.h"
#include "NRUTIL.h"


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
	    print_gmsh = (char*) malloc(3 * sizeof(char *));
	    fname      = (char*) malloc(23 * sizeof(char));
	    input_inp  = (char*) malloc(8 * sizeof(char));
  
  	    
	    outfile_msh_vtk  = (char*) malloc(96 * sizeof(char *));
	    outfile_msh_gmsh = (char*) malloc(96 * sizeof(char *));
	    outfile_msh_alya = (char*) malloc(96 * sizeof(char *));
  		external_grid_file_name = (char*) malloc(255 * sizeof(char *));
  		*external_grid_file_name = 0;
	    
	    parameters = dvector(0,NUMBER_OF_PARAMETERS);

	    //Allocation of the strings of the above lines:
	    for(i=0; i<PROB_ENTRIES; i++)
			problem[i] = (char *) malloc(32 * sizeof(char *));
      
	    printf(" # Memory allocated (flag 0): various input-related strings\n");
	}
    if(flag == 1)
	{
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

	    printf(" # Memory allocated (flag 1):\n");
	    printf(" # \t\t COORDS[nnodes+1][nsd+1]\n");
	    printf(" # \t\t COORDS1d[nnodes*(nsd+1)]\n");
	    printf(" # \t\t CONN[nelem+1][EL_NODES+1]\n");
	    printf(" # \t\t ELTYPE[nelem+1]\n");
	    printf(" # \t\t BDYFLAG[nnodes+1][nsd+1]\n");
	    
	}

    if(flag == 2)
	{     
	    //Allocate space for the array of boundary coordinates to be read in input:
	    BDY_COORDS = dmatrix(0,NBDY_NODES,0,nsd);	   
	    printf(" # Memory allocated (flag 2): BDY_COORDS[NBDY_NODES+1][nsd+1] %d\n", nsd);
	}

    if(flag == 3)
	{
	    rarray = dvector(0,nnodes);	   
	    printf(" # Memory allocated (flag 3): rarray[nnodes+1]\n");
	}

    if(flag == 4)
	{     
	    //Allocate space for the 3d matrices x,y,z
	    x = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);
	    y = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);
	    z = d3tensor(0,nnodesx,0,nnodesy,0,nnodesz);

	    printf(" # Memory allocated (flag 4):\n");
	    printf(" # \t\t x[nnodesx+1][nnodesy+1][nnodesz+1]\n");
	    printf(" # \t\t y[nnodesx+1][nnodesy+1][nnodesz+1]\n");
	    printf(" # \t\t z[nnodesx+1][nnodesy+1][nnodesz+1]\n");
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
  
  		free(external_grid_file_name);
	    free(outfile_msh_vtk);
	    free(outfile_msh_alya);
  
	    for(i=0; i<PROB_ENTRIES; i++)
		free(problem[i]);

	    printf(" # Freed memory (flag 0) for input-related strings.\n");
	}
    if(flag == 1)
	{
      
	    free_dvector(parameters,0,NUMBER_OF_PARAMETERS);
	    free_dvector(COORDS1d,  0,nnodes*(nsd+1));
	    free_ivector(ELTYPE,    0,nelem+1);
	    free_dmatrix(COORDS,    0,nnodes,     0,nsd+1);
	    free_imatrix(CONN,      0,nelem,      0,EL_NODES+1);
	    // free_dmatrix(BDY_COORDS,0,NBDY_NODES, 0,nsd);

	    printf(" # Freed memory (flag 1):\n");
	    printf(" # \t\t free(COORDS)\n");
	    printf(" # \t\t free(COORDS1d)\n");
	    printf(" # \t\t free(CONN)\n");
	    printf(" # \t\t free(ELTYPE)\n");
	    printf(" # \t\t free(BDYFLAG)\n");
	}
    if(flag == 2)
	{     
	    free_dmatrix(BDY_COORDS, 0,NBDY_NODES, 0,nsd);
	    printf(" # Freed memory (flag 2): free(BDY_COORDS)\n");
	}

    if(flag == 3)
	{
	    free_dvector(rarray, 0,nnodesx*nnodesy);
	    printf(" # Freed memory (flag 3): free(rarray)\n");
	}
    if(flag == 4)
	{     
	    //Allocate space for the array of boundary coordinates to be read in input:
	    free_d3tensor(x,0,nnodesx,0,nnodesy,0,nnodesz);
	    free_d3tensor(y,0,nnodesx,0,nnodesy,0,nnodesz);
	    free_d3tensor(z,0,nnodesx,0,nnodesy,0,nnodesz);

	    printf(" # Freed memory (flag 4):\n");
	    printf(" # \t\t free(x)\n");
	    printf(" # \t\t free(y)\n");
	    printf(" # \t\t free(z)\n");
	    
	}

    return 0;    
}

