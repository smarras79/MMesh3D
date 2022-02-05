/*
 * This function writes the outputs to file (vtk, alya, etc.):
 *
 * Simone Marras, Aug 2012
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "GLOBAL_VARS.h"
#include "MYDEFINE.h"
#include "PRINT.h"
#include "WRITE_OUTPUT.h"

int WRITE_OUTPUT(int rank)
{
	
    //Write scalars to strings:
    sprintf(nodesx, "%d", nnodesx);
    sprintf(nodesy, "%d", nnodesy);
    sprintf(nodesz, "%d", nnodesz); 
    sprintf(nodes,  "%d", nnodes);
    sprintf(elorder,"%d", nop);
    sprintf(nel,"%d", nelem);
    sprintf(mpiprocess,"%d", rank);
    
  
    /*****************************************************************************************
     * Write the VTK plot file:
     * To do so we need to re-calculate the total nnodes and nelem of the global final grid
     *****************************************************************************************/
    //VTK output file name:
    if( !strcmp(print_vtk, "yes") || !strcmp(print_vtk,"y") || !strcmp(print_vtk,"vtk"))
	{
	    strcpy(outfile_msh_vtk,"VTK_mesh_MPIRANK_");
	    strcat(outfile_msh_vtk,mpiprocess);
	    strcat(outfile_msh_vtk,"_npoints_");
	    strcat(outfile_msh_vtk,nodes);
	    strcat(outfile_msh_vtk,"_nelem_");
	    strcat(outfile_msh_vtk,nel);
	    strcat(outfile_msh_vtk,"_nop_");
	    strcat(outfile_msh_vtk,elorder);
	    strcat(outfile_msh_vtk, ".vtk");
	  
	    dPRINT_UNSTRUCT_GRID_VTK(outfile_msh_vtk,  1, 1);
      
	    if(rank == 0){		
		printf(" #\n");
		printf(" # The following file was created:\n");
		printf(" # %s:\n", outfile_msh_vtk);
	    }
	}
  
    /*****************************************************************************************
     * Write the GMSH/ABAQUS plot file:
     *****************************************************************************************/
  
    //GMSH/ABAQUS output file name:
    if( !strcmp(print_gmsh, "yes") || !strcmp(print_gmsh,"y") || !strcmp(print_gmsh,"gmsh"))
	{
	    strcpy(outfile_msh_gmsh,"GMSH_mesh_MPIRANK_");
	    strcat(outfile_msh_gmsh,mpiprocess);
	    strcat(outfile_msh_gmsh,"_");
	    strcat(outfile_msh_gmsh,nodesx);
	    strcat(outfile_msh_gmsh,"x");
	    strcat(outfile_msh_gmsh,nodesy);
	    strcat(outfile_msh_gmsh,"x");
	    strcat(outfile_msh_gmsh,nodesz);
	    strcat(outfile_msh_gmsh,"x");
	    strcat(outfile_msh_gmsh,elorder);
	    strcat(outfile_msh_gmsh, ".gmsh");
      
	    dPRINT_UNSTRUCT_GRID_GMSH(outfile_msh_gmsh,  1, 1);

	    if(rank == 0){		
		printf(" #\n");
		printf(" # The following file was created:\n");
		printf(" # %s:\n", outfile_msh_gmsh);
	    }
	}


    /*
     * WRITE THE CONNECTIVITY FILE TO BE USED BY METIS (m2gmetis) 
     * that transforms the grid format into the graph that metis needs
     */
    strcpy(outfile_msh_gmsh,"CONN_mesh_MPIRANK_");
    strcat(outfile_msh_gmsh,mpiprocess);
    strcat(outfile_msh_gmsh,"_nelem_");
    strcat(outfile_msh_gmsh,nel);
    strcat(outfile_msh_gmsh,"_nop_");
    strcat(outfile_msh_gmsh,elorder);
    strcat(outfile_msh_gmsh, ".dat");
  
    dPRINT_UNSTRUCT_GRID_CONN(outfile_msh_gmsh,  1, 1);

    if(rank == 0){		
	printf(" #\n");
	printf(" # The following file was created:\n");
	printf(" # %s:\n", outfile_msh_gmsh);
    }

    return 0;
}
