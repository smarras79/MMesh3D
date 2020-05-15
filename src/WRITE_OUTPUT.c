/*
 * This function writes the outputs to file (vtk, alya, etc.):
 *
 * Simone Marras, Aug 2012
 */
#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

int WRITE_OUTPUT(int rank)
{
  //Call MPI rank:
  

  //Write scalars to strings:
  sprintf(nodesx, "%d", nnodesx);
  sprintf(nodesy, "%d", nnodesy);
  sprintf(nodesz, "%d", nnodesz); 
  sprintf(elorder,"%d", nop);
  sprintf(mpiprocess,"%d", rank);
  
  if( !strcmp(print_alya, "yes") || !strcmp(print_alya,"y") || !strcmp(print_alya,"alya"))
    {
      
      //Mesh ALYA output file:
      strcpy(outfile_msh_alya, "3D_mesh_");
      strcat(outfile_msh_alya,problem[2]);
      strcat(outfile_msh_alya,"_");
      strcat(outfile_msh_alya,mpiprocess);
      strcat(outfile_msh_alya,"_");
      strcat(outfile_msh_alya,nodesx);
      strcat(outfile_msh_alya,"x");
      strcat(outfile_msh_alya,nodesy);
      strcat(outfile_msh_alya,"x");
      strcat(outfile_msh_alya,nodesz);
      strcat(outfile_msh_alya,"x");
      strcat(outfile_msh_alya,elorder);
      strcat(outfile_msh_alya, ".dom.dat");

      dGRID2ALYAFILE(outfile_msh_alya, COORDS, CONN, nnodes, nelem, 0);

      if(rank == 0){
	printf(" # The following file was created:\n");
	printf(" # %s:\n\n", outfile_msh_alya);
      }
      
    }
  
  /*****************************************************************************************
   * Write the VTK plot file:
   * To do so we need to re-calculate the total nnodes and nelem of the global final grid
   *****************************************************************************************/
  
  //VTK output file name:
  if( !strcmp(print_vtk, "yes") || !strcmp(print_vtk,"y") || !strcmp(print_vtk,"vtk"))
    {
      strcpy(outfile_msh_vtk,"VTK_mesh_");
      strcat(outfile_msh_vtk,problem[2]);
      strcat(outfile_msh_vtk,"_");
      strcat(outfile_msh_vtk,mpiprocess);
      strcat(outfile_msh_vtk,"_");
      strcat(outfile_msh_vtk,nodesx);
      strcat(outfile_msh_vtk,"x");
      strcat(outfile_msh_vtk,nodesy);
      strcat(outfile_msh_vtk,"x");
      strcat(outfile_msh_vtk,nodesz);
      strcat(outfile_msh_vtk,"x");
      strcat(outfile_msh_vtk,elorder);
      strcat(outfile_msh_vtk, ".vtk");
      
      dPRINT_UNSTRUCT_GRID_VTK(outfile_msh_vtk,  1, 1);
      
      if(rank == 0){
	printf(" # The following file was created:\n");
	printf(" # %s:\n\n", outfile_msh_vtk);
      }

    }
  
  /*****************************************************************************************
   * Write the GMSH/ABAQUS plot file:
   *****************************************************************************************/
  
  //GMSH/ABAQUS output file name:
  if( !strcmp(print_gmsh, "yes") || !strcmp(print_gmsh,"y") || !strcmp(print_gmsh,"gmsh"))
    {
      strcpy(outfile_msh_gmsh,"GMSH_mesh_");
      strcat(outfile_msh_gmsh,problem[2]);
      strcat(outfile_msh_gmsh,"_");
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
	printf(" # The following file was created:\n");
	printf(" # %s:\n\n", outfile_msh_gmsh);
      }
    }


  /*
   * WRITE THE CONNECTIVITY FILE TO BE USED BY METIS (m2gmetis) 
   * that transforms the grid format into the graph that metis needs
   */
  strcpy(outfile_msh_gmsh,"CONN_mesh_");
  strcat(outfile_msh_gmsh,problem[2]);
  strcat(outfile_msh_gmsh,"_");
  strcat(outfile_msh_gmsh,mpiprocess);
  strcat(outfile_msh_gmsh,"_");
  strcat(outfile_msh_gmsh,nodesx);
  strcat(outfile_msh_gmsh,"x");
  strcat(outfile_msh_gmsh,nodesy);
  strcat(outfile_msh_gmsh,"x");
  strcat(outfile_msh_gmsh,nodesz);
  strcat(outfile_msh_gmsh,"x");
  strcat(outfile_msh_gmsh,elorder);
  strcat(outfile_msh_gmsh, ".dat");
  
  dPRINT_UNSTRUCT_GRID_CONN(outfile_msh_gmsh,  1, 1);

  if(rank == 0){
    printf(" # The following file was created:\n");
    printf(" # %s:\n\n", outfile_msh_gmsh);
  }

  return 0;
}
