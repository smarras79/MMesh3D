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
    if(flag == 5)
	{
	    CONN_BDY_FACES = imatrix(0,nbdy_faces, 0,(nop+1)*(nop+1)); //For hexa only for now
	    for(int ibdy_face=0; ibdy_face<=nbdy_faces; ibdy_face++){
		for(inode=0; inode<((nop+1)*(nop+1)); inode++){
		    CONN_BDY_FACES[ibdy_face][inode] = 0;
		}
	    }
	    
	    printf(" # Memory allocated (flag 5): CONN_BDY_EDGES[nbdy_faces][(nop+1)*(nop+1)]\n");
	}
    if(flag == 6)
	{
	    FACE_LtoG    = imatrix (0,nelem-1, 0, 5);
	    FACE_in_ELEM = i3tensor(0,nelem-1, 0, 5, 0,1);
	    conn_edge_el = i3tensor(0,nelem-1, 0,11, 0,1);
	    conn_face_el = i3tensor(0,nelem-1, 0, 5, 0,3);
	    
	    printf(" # Memory allocated (flag 6):\n");
	    printf(" # \t\t FACE_LtoG[nelem][6]\n");
	    printf(" # \t\t FACE_in_ELEM[nelem][12][2]\n");
	    printf(" # \t\t conn_edge_el[nelem][12][2]\n");
	    printf(" # \t\t conn_face_el[nelem][6][4]\n");
	}
    if(flag == 7)
	{
	    CONN_FACE = imatrix(0,nfaces, 0, (nop+1)*(nop+1)); //For hexa only for now
	    for(int iface=0; iface<nfaces; iface++){
		for(int inode=0; inode<((nop+1)*(nop+1)); inode++){
		    CONN_FACE[iface][inode] = 0;
		}
	    }
	    printf(" # Memory allocated (flag 7): CONN_FACES[%d][%d]\n", nfaces, (nop+1)*(nop+1));
	}
    if(flag == 8)
	{
	    CONN_EDGE = imatrix(0,nedges, 0, (nop+1)); //For hexa only for now
	    for(int iedge=0; iedge<nedges; iedge++){
		for(int inode=0; inode<(nop+1); inode++){
		    CONN_EDGE[iedge][inode] = 0;
		}
	    }
	    printf(" # Memory allocated (flag 8): CONN_EDGES[%d][%d]\n", nedges, (nop+1));
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
    if(flag == 5)
	{
	    free_imatrix(CONN_BDY_FACES, 0,nbdy_faces, 0,((nop+1)*(nop+1))); //For hexa only for now
	    printf(" # Freed memory (flag 5): free(CONN_BDY_FACE)\n");
	}
    if(flag == 6)
	{
	    free_imatrix(FACE_LtoG, 0,nelem-1, 0, 5);
	    free_i3tensor(FACE_in_ELEM,  0,nelem-1, 0, 5, 0,1);
	    free_i3tensor(conn_edge_el,  0,nelem-1, 0,11, 0,1);
	    free_i3tensor( conn_face_el, 0,nelem-1, 0, 5, 0,3);
	     
	    printf(" # Freed memory (flag 6):\n");
	    printf(" # \t\t free(FACE_LtoG)\n");
	    printf(" # \t\t free(FACE_in_ELEM)\n");
	    printf(" # \t\t free(conn_edge_el)\n");
	    printf(" # \t\t free(conn_face_el)\n");
	}
     if(flag == 7)
	{
	    free_imatrix(CONN_FACE, 0,nfaces, 0,((nop+1)*(nop+1))); //For hexa only for now
	    printf(" # Freed memory (flag 7): free(CONN_FACE)\n");
	}
     if(flag == 8)
	{
	    free_imatrix(CONN_EDGE, 0,nedges, 0, (nop+1)); //For hexa only for now
	    printf(" # Freed memory (flag 8): free(CONN_EDGES)\n");
	}
    
     
    return 0;    
}

