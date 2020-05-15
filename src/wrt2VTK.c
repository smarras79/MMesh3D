/*
 * This function writes a VTK file for data on
 * an unstructured grid.
 *
 * simone.marras@gmail.com
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void wrt2VTK(double **COORDS, int **CONN, int *ELTYPE, int nelem, int nnodes,
	     char *file_name, int nsd)
{  
  int i,ie;
  int nsize, nfield_data;
  
  FILE *file_id;
  
  nsize = 5*nelem;
  
  file_id = fopen(file_name, "w");
  if(file_id == NULL)
    printf("The file %s could not be open\n", file_name);
  else{
    
    //Open VTK file
    fprintf(file_id, "# vtk DataFile Version 2.0\n");
    fprintf(file_id, "INIT data\n");
    fprintf(file_id, "ASCII\n");
    fprintf(file_id, "DATASET UNSTRUCTURED_GRID\n");
    
    //Write coordinates
    if( nsd == 2 )
      {
	fprintf(file_id, "POINTS %d double\n", nnodes);
	for(i=1; i<=nnodes; i++)
	  {
	    fprintf(file_id, " %lf %lf %lf\n", COORDS[i][1], COORDS[i][2], 0.0);
	  }
	
	//Write coonectivity
	fprintf(file_id, "CELLS %d %d\n", nelem,nsize);
	for(ie=1; ie<=nelem; ie++)
	  {
	    fprintf(file_id, " %d %d %d %d %d\n", ELTYPE[ie], CONN[ie][2]-1, CONN[ie][3]-1, CONN[ie][4]-1, CONN[ie][5]-1);
	  }
	
	fprintf(file_id, "CELL_TYPES %d\n", nelem);
	for(ie=1; ie<=nelem; ie++)
	  fprintf(file_id, " %d\n", 9);

      }
    else if( nsd == 3 )
      {
	fprintf(file_id, "POINTS %d double\n", nnodes);
	for(i=1; i<=nnodes; i++)
	  {
	    fprintf(file_id, " %lf %lf %lf\n", COORDS[i][1], COORDS[i][2], COORDS[i][3]);
	  }
	
	//Write coonectivity
	fprintf(file_id, "CELLS %d %d\n", nelem,nsize);
	for(ie=1; ie<=nelem; ie++)
	  {
	    fprintf(file_id, " %d %d %d %d %d\n", ELTYPE[ie], CONN[ie][2]-1, CONN[ie][3]-1, CONN[ie][4]-1, CONN[ie][5]-1);
	  }
	
	fprintf(file_id, "CELL_TYPES %d\n", nelem);
	for(ie=1; ie<=nelem; ie++)
	  {
	    if( ELTYPE[ie] == 8 ) //HEXA
	      fprintf(file_id, " %d\n", 12);
	    else if ( ELTYPE[ie] == 6 ) //WEDGES
	      fprintf(file_id, " %d\n", 13);
	  }
    
      }
    fclose(file_id);
  }
  return;
}

