/*PRINT.c*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "visit_writer.h"
#include "PRINT.h"
#include "MYDEFINE.h"
#include "GLOBAL_VARS.h"

void MAT2F_f2DMAT(char *mat_name, float **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *mat_id;
	
    if((mat_id = fopen(mat_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    for(j=min_cols; j<=max_cols; j++)
		fprintf(mat_id, " %f", A[i][j]);
	    fprintf(mat_id, "\n");
	}
	//printf("\n");
    }
    fclose(mat_id);
		
    return;
}

void MAT2F_d2DMAT(char *mat_name, double **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *mat_id;
	
    if((mat_id = fopen(mat_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    for(j=min_cols; j<=max_cols; j++)
		fprintf(mat_id, " %f", A[i][j]);
	    fprintf(mat_id, "\n");
	}
	//printf("\n");
    }
    fclose(mat_id);
		
    return;
}

void MAT2F_i2DMAT(char *mat_name, int **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *mat_id;
	
    if((mat_id = fopen(mat_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    for(j=min_cols; j<=max_cols; j++)
		fprintf(mat_id, " %i", A[i][j]);
	    fprintf(mat_id, "\n");
	}
	//printf("\n");
    }
    fclose(mat_id);
		
    return;
}

void VECT2F_fVECT(char *vect_name, float *V, int min_rows, int max_rows)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *vect_id;
	
    if((vect_id = fopen(vect_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    fprintf(vect_id, " %f\n", V[i]);
	}
    }
    fclose(vect_id);
	
    return;
}

void VECT2F_dVECT(char *vect_name, double *V, int min_rows, int max_rows)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *vect_id;
	
    if((vect_id = fopen(vect_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    fprintf(vect_id, " %f\n", V[i]);
	}
    }
    fclose(vect_id);
		
    return;
}

void VECT2F_iVECT(char *vect_name, int *V, int min_rows, int max_rows)
{
    //Function to WRITE a matrix of doubles to a text file
    int i,j;
    FILE *vect_id;
	
    if((vect_id = fopen(vect_name, "w")) == NULL)
	printf("The file could not be open\n");
    else{
	for(i=min_rows; i<=max_rows; i++){
	    fprintf(vect_id, " %i\n", V[i]);
	}
    }
    fclose(vect_id);
		
    return;
}

void VIEW_f2DMAT(char *mat_name, float **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to view a matrix of doubles
    int i,j;
	
    printf("%s \n", mat_name);
    for(i=min_rows; i<=max_rows; i++){
	for(j=min_cols; j<=max_cols; j++)
	    printf(" %.2f", A[i][j]);
	printf("\n");
    }
    printf("\n");
		
    return;
}

void VIEW_d2DMAT(char *mat_name, double **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to view a matrix of doubles
    int i,j;
	
    printf("%s \n", mat_name);
    for(i=min_rows; i<=max_rows; i++){
	for(j=min_cols; j<=max_cols; j++)
	    printf(" %f", A[i][j]);
	printf("\n");
    }
    printf("\n");
		
    return;
}

void VIEW_i2DMAT(char *mat_name, int **A, int min_rows, int max_rows, int min_cols, int max_cols)
{
    //Function to view a matrix of integers
    int i,j;
	
    printf("%s \n", mat_name);
    for(i=min_rows; i<=max_rows; i++){
	for(j=min_cols; j<=max_cols; j++)
	    printf(" %d", A[i][j]);
	printf("\n");
    }
    printf("\n");
		
    return;
}

void VIEW_i3DMAT(char *mat_name, int ***A, int min_rows, int max_rows, int min_cols, int max_cols, int min_depth, int max_depth)
{
    //Function to view a matrix of integers
    int i,j;
	
    printf("%s \n", mat_name);    
    for(int i=min_rows; i<=max_rows; i++){
	for(int j=min_cols; j<=max_cols; j++){
	    for(int k=min_depth; k<=max_depth; k++){
		printf(" %d", A[i][j][k]);
	    }
	    printf("\n");
	}
	printf("\n");
    }
    
    return;
}

void VIEW_fVECT(char *vect_name, float *V, int min_rows, int max_rows)
{
    //Function to view a 1D array of doubles
    int i;
			
    printf("%s \n", vect_name);
    for(i=min_rows; i<=max_rows; i++)
	printf(" %f\n", V[i]);
    printf("\n");

    return;
}


void VIEW_dVECT(char *vect_name, double *V, int min_rows, int max_rows)
{
    //Function to view a 1D array of doubles
    int i;
			
    printf("%s \n", vect_name);
    for(i=min_rows; i<=max_rows; i++)
	printf(" %f\n", V[i]);
    printf("\n");

    return;
}

void VIEW_iVECT(char *vect_name, int *V, int min_rows, int max_rows)
{
    //Function to view a 1D array of integers
    int i;
	
    printf("%s \n", vect_name);
    for(i=min_rows; i<=max_rows; i++)
	printf(" %d\n", V[i]);
    printf("\n");
			
    return;
}

void PRINT(float **u, float **v, float **P, int imax, int jmax)
{

    int i,j;

    //PRINT U,V,P
    printf("U[i][j]: \n");
    for(i=0; i<=imax; i++){
	for(j=0; j<=jmax+1; j++)
	    printf("%f ", u[i][j]);
			
	printf(" \n");
    }
	
    printf("V[i][j]: \n");
    for(i=0; i<=imax+1; i++){
	for(j=0; j<=jmax; j++)
	    printf("%f ", v[i][j]);
			
	printf(" \n");
    }
		
    printf("P[i][j]: \n");
    for(i=0; i<=imax+1; i++){
	for(j=0; j<=jmax+1; j++)
	    printf("%f ", P[i][j]);
			
	printf(" \n");
    }
}

/****************************************/
/**PRINT_VECT_VTK.c**********************/
/****************************************/
void PRINT_VECT_VTK(int imax, int jmax, int kmax, char *vector_field_name, float **U, float **V, float **W, float *x, float *y, float *z)
{
    int i,j,k;
    int nx,ny,nz;
	
    FILE *vect_field;
	
    nx = imax+1;
    ny = jmax+1;
    nz = 1;
	
    /* File grid.vtk*/
    if((vect_field = fopen(vector_field_name, "w")) == NULL)
	printf("The grid file could not be open\n");
    else{
	fprintf(vect_field, "# vtk DataFile Version 2.0\n");
	fprintf(vect_field, "Sample rectilinear grid\n");
	fprintf(vect_field, "ASCII\n");
	fprintf(vect_field, "DATASET RECTILINEAR_GRID\n");
	fprintf(vect_field, "DIMENSIONS %d %d %d\n", ny, nx, nz);
	fprintf(vect_field, "Y_COORDINATES %d float\n", ny);
	for(j=0; j<=ny-1; j++)
	    fprintf(vect_field, "%f\n", y[j]);
	fprintf(vect_field, "X_COORDINATES %d float\n", nx);	
	for(i=0; i<=nx-1; i++)
	    fprintf(vect_field, "%f\n", x[i]);	
	//fprintf(vect_field, "Z_COORDINATES %d float\n", nz);
	//fprintf(vect_field, "%f\n", 0.0*z);
	fprintf(vect_field, "POINT_DATA %d\n", (nx)*(ny)*(nz));
	fprintf(vect_field, "VECTORS vectors float\n");
	for(j=0; j<=ny-1; j++)
	    for(i=0; i<=nx-1; i++){
		fprintf(vect_field, "%f\t %f\t %f\n", U[i][j], V[i][j], 0.0);	
	    }
    }
    fclose(vect_field);

    return;
}


/****************************************/
/**PRINT_GRID_VTK.c**********************/
/****************************************/
void PRINT_GRID_VTK(char *grid_name, int imax, int jmax, int kmax, float *x, float *y, float *z)
{
    int i,j,k;
    int nx,ny,nz;
	
    FILE *grid;
	
    nx = imax+1;
    ny = jmax+1;
    nz = kmax+1;
	
    //Control if the file of the grid must be written or not
    //and checked before going on with the computations
	
    /* File grid.vtk*/
    if((grid = fopen(grid_name, "w")) == NULL)
	printf("The grid file could not be open\n");
    else{
	fprintf(grid, "# vtk DataFile Version 2.0\n");
	fprintf(grid, "Sample rectilinear grid\n");
	fprintf(grid, "ASCII\n");
	fprintf(grid, "DATASET RECTILINEAR_GRID\n");
	fprintf(grid, "DIMENSIONS %d %d %d\n", ny, nx, nz);
	fprintf(grid, "Z_COORDINATES %d float\n", nz);
	for(k=0; k<=nz-1; k++)
	    fprintf(grid, "%f\n", z[k]);
	fprintf(grid, "Y_COORDINATES %d float\n", ny);
	for(j=0; j<=ny-1; j++)
	    fprintf(grid, "%f\n", y[j]);
	fprintf(grid, "X_COORDINATES %d float\n", nx);	
	for(i=0; i<=nx-1; i++)
	    fprintf(grid, "%f\n", x[i]);
	fprintf(grid, "POINT_DATA %d\n", (nx)*(ny)*(nz));
    }
    fclose(grid);
	
    // printf("\nWARNING!!!\nYou chose 'YES' in the input '*print_grid' to view the grid file\n");
    // printf("Follow the instruction to view it and go back to the computations\n");
    // printf("\n 1) Open PARAVIEW from X11 and Open the grid.vtk\n");
    // printf(" 2) Check if the generated grid is as expected\n");
    // printf(" 3) If yes, change the input *print_grid to 'NO' and RECOMPILE and execute the program\n\n");

    return;
}

/****************************************/
/**PRINT_STRUCT_GRID_VTK.c***************/
/****************************************/
void PRINT_STRUCT_GRID_VTK(char *grid_name, int imax, int jmax, int kmax, float *x, float *y, float *z)
{
    int i,j,k,id=0;
    int nx,ny,nz, ntot;
	
    FILE *grid;
	
    nx = imax+1;
    ny = jmax+1;
    nz = kmax+1;
    ntot = nx*ny*nz;
	
    //Control if the file of the grid must be written or not
    //and checked before going on with the computations
	
    /* File grid.vtk*/
    if((grid = fopen(grid_name, "w")) == NULL)
	printf("The grid file could not be open\n");
    else{
	fprintf(grid, "# vtk DataFile Version 2.0\n");
	fprintf(grid, "myMeshGen VTK file\n");
	fprintf(grid, "ASCII\n");
	fprintf(grid, "DATASET STRUCTURED_GRID\n");
	fprintf(grid, "DIMENSIONS %d %d %d\n", ny, nx, nz);
	fprintf(grid, "POINTS %d float\n", ntot);
	for(k=0; k<=nz-1; k++){
	    for(j=0; j<=ny-1; j++)
		for(i=0; i<=nx-1; i++)
		    fprintf(grid, "%f %f %f\n", x[i], y[j], z[k]);
	}
	fclose(grid);
    }
    // printf("\nWARNING!!!\nYou chose 'YES' in the input '*print_grid' to view the grid file\n");
    // printf("Follow the instruction to view it and go back to the computations\n");
    // printf("\n 1) Open PARAVIEW from X11 and Open the grid.vtk\n");
    // printf(" 2) Check if the generated grid is as expected\n");
    // printf(" 3) If yes, change the input *print_grid to 'NO' and RECOMPILE and execute the program\n\n");

    return;
}

/***************************************************************
 * PRINT_UNSTRUCT_GRID_VTK.c
 * this function writes the unstructured grid into a VTK file
 * It uses the function write_unstructured_mesh_only from
 * the "modified" library visit_writer that comes with this package.
 *
 *
 * See the visit_writer.h for the copyright issues of the library.
 * simone.marras@gmail.com
 ****************************************************************/
void dPRINT_UNSTRUCT_GRID_VTK(char *grid_name, int array_numb, int coords_ncolumns)

{
    unsigned int i, j, k;
    	
    //visit writer wants the COORDS and CONN elements in a sequential 1D array
    //int elem_visit_types[nelem+1];

    //The way I wrote this loops may seems more
    //complex than what it really is; but it is poweful
    //because it can take in input any type of connectivity in terms of max_elnodes;
    //If you do it by hand, you can figure what is doing!
    for( i = array_numb; i<=nelem-1+array_numb; i++){
	//printf("nelem = %d\n", i);
	if(ELTYPE[i] == 3){
	    ELTYPE[i-1] = VISIT_TRIANGLE;}
	else if(ELTYPE[i] == 4){
	    ELTYPE[i-1] = VISIT_QUAD;}
	else if(ELTYPE[i] == 8){
	    ELTYPE[i-1] = VISIT_HEXAHEDRON;}
	else if(ELTYPE[i] == 6){
	    ELTYPE[i-1] = VISIT_WEDGE;
	}
    }
    
    /* Pass the mesh and data to visit_writer. */    
    dwrt2VTK(grid_name);
        
    return;
}

//As above but in double precision:
void dPRINT_UNSTRUCT_GRID_GMSH(char *grid_name, int array_numb, int coords_ncolumns)
{
 
    unsigned int i, j, k;
	
    //visit writer wants the COORDS and CONN elements in a sequential 1D array
    //int elem_visit_types[nelem+1];

    //The way I wrote this loops may seems more
    //complex than what it really is; but it is poweful
    //because it can take in input any type of connectivity in terms of max_elnodes;
    //If you do it by hand, you can figure what is doing!
    for( i = array_numb; i<=nelem-1+array_numb; i++){
	//printf("nelem = %d\n", i);
	if(ELTYPE[i] == 3){
	    ELTYPE[i-1] = VISIT_TRIANGLE;}
	else if(ELTYPE[i] == 4){
	    ELTYPE[i-1] = VISIT_QUAD;}
	else if(ELTYPE[i] == 8){
	    ELTYPE[i-1] = VISIT_HEXAHEDRON;}
	else if(ELTYPE[i] == 6){
	    ELTYPE[i-1] = VISIT_WEDGE;
	}
    }
	
    /* Pass the mesh and data to visit_writer. */
    dwrt2GMSH(grid_name);

    return;
}

//As above but in double precision:
void dPRINT_UNSTRUCT_GRID_CONN(char *grid_name, int array_numb, int coords_ncolumns)
{
 
    unsigned int i, j, k;
	
    //visit writer wants the COORDS and CONN elements in a sequential 1D array
    //int elem_visit_types[nelem+1];

    //The way I wrote this loops may seems more
    //complex than what it really is; but it is poweful
    //because it can take in input any type of connectivity in terms of max_elnodes;
    //If you do it by hand, you can figure what is doing!
    for( i = array_numb; i<=nelem-1+array_numb; i++){
	//printf("nelem = %d\n", i);
	if(ELTYPE[i] == 3){
	    ELTYPE[i-1] = VISIT_TRIANGLE;}
	else if(ELTYPE[i] == 4){
	    ELTYPE[i-1] = VISIT_QUAD;}
	else if(ELTYPE[i] == 8){
	    ELTYPE[i-1] = VISIT_HEXAHEDRON;}
	else if(ELTYPE[i] == 6){
	    ELTYPE[i-1] = VISIT_WEDGE;
	}
    }
  
    /* Pass the mesh and data to visit_writer. */
    dwrt2CONN(grid_name);

    return;
}



void dwrt2VTK_nx_ny_nz(char *file_name)
{  
    int i,j,ie,inode;
    int u,v,w;
    int elu,elv,elw, iel;
    int nsize, nfield_data;
  
    FILE *file_id;
  
  
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
		fprintf(file_id, "POINTS %d float\n", nnodes);
		for(i=1; i<=nnodes; i++)
		    {
			fprintf(file_id, " %.12f %.12f %.12f\n", COORDS[i][2], COORDS[i][3], 0.0);
		    }
	
		//Write coonectivity
		nsize = 5*nelem;
		fprintf(file_id, "CELLS %d %d\n", nelem, nsize);
		for(ie=1; ie<=nelem; ie++)
		    {
			if( ELTYPE[ie] == VISIT_QUAD )
			    fprintf(file_id, " %d %d %d %d %d\n", ELTYPE[ie], MAPL2G[ie][2]-1, MAPL2G[ie][3]-1, MAPL2G[ie][4]-1, MAPL2G[ie][5]-1);
			else if( ELTYPE[ie] == VISIT_TRIANGLE )
			    fprintf(file_id, " %d %d %d %d\n", ELTYPE[ie], MAPL2G[ie][2]-1, MAPL2G[ie][3]-1, MAPL2G[ie][4]-1);
		    }
	
		fprintf(file_id, "CELL_TYPES %d\n", nelem);
		for(ie=1; ie<=nelem; ie++)
		    fprintf(file_id, " %d\n", ELTYPE[ie]);
	
	    }
	else if( nsd == 3 )
	    {
		fprintf(file_id, "POINTS %d double\n", nnodes);

		i=1;
		j=1;
		for (w=1; w<=nnodesz; w++){
		    for (v=1; v<=nnodesy; v++){
			for (u=1; u<=nnodesx; u++){
	    
			    fprintf(file_id, " %.12f %.12f %.12f\n", COORDS[i][2], COORDS[i][3], COORDS[i][4]);
	      
			    i++;
			}
		    }
		}
	
		//COORDS --> COORDS1d:
		i=1;
		j=1;
		for (w=1; w<=nnodesz; w++){
		    for (v=1; v<=nnodesy; v++){
			for (u=1; u<=nnodesx; u++){
	      
			    COORDS1d[j]   = COORDS[i][2];
			    COORDS1d[j+1] = COORDS[i][3];
			    COORDS1d[j+2] = COORDS[i][4];
	      
			    i++;
			    j=j+3;
			}
		    }
		}
	
		/* Write the Coordinates 1D array to file */
		VECT2F_dVECT("COORDS1d.dat", COORDS1d, 1,nnodes*3);
	
		//Write coonectivity
		if( ELTYPE[1] == VISIT_HEXAHEDRON )
		    nsize = 9*nelem;
		else if( ELTYPE[1] == VISIT_WEDGE )
		    nsize = 7*nelem;

		fprintf(file_id, "CELLS %d %d\n", nelem, nsize);
		for(ie=1; ie<=nelem; ie++)
		    {
			if( ELTYPE[ie-1] == VISIT_HEXAHEDRON )
			    fprintf(file_id, " %d %d %d %d %d %d %d %d %d\n", 8, MAPL2G[ie][2]-1, MAPL2G[ie][3]-1, MAPL2G[ie][4]-1, MAPL2G[ie][5]-1, MAPL2G[ie][6]-1, MAPL2G[ie][7]-1, MAPL2G[ie][8]-1, MAPL2G[ie][9]-1);
			else if( ELTYPE[ie-1] == VISIT_WEDGE )
			    fprintf(file_id, " %d %d %d %d %d %d %d\n", 6, MAPL2G[ie][2]-1, MAPL2G[ie][3]-1, MAPL2G[ie][4]-1, MAPL2G[ie][5]-1, MAPL2G[ie][6]-1, MAPL2G[ie][7]-1);
		    }
	
		fprintf(file_id, "CELL_TYPES %d\n", nelem);
		for(ie=0; ie<nelem; ie++)
		    fprintf(file_id, " %d\n", VISIT_HEXAHEDRON);

	    }
    }
    fclose(file_id);
  
    return;
}


void dwrt2VTK(char *file_name)
{  
    int i, j, ie, inode;
    int ii, jj, kk;
    int u, v, w;
    int elu, elv, elw, iel;
    int nsize, nfield_data;
    int ncells;
    
    FILE *file_id;
    
    file_id = fopen(file_name, "w");
    if(file_id == NULL)
	printf("The file %s could not be open\n", file_name);
    else{
    
	//Open VTK file
	fprintf(file_id, "# vtk DataFile Version 2.0\n");
	fprintf(file_id, "MESH data\n");
	fprintf(file_id, "ASCII\n");
	fprintf(file_id, "DATASET UNSTRUCTURED_GRID\n");

	if (nop < 2) {
	    //Write coordinates
	    if( nsd == 2 ) {
		fprintf(file_id, "POINTS %d float\n", nnodes);
		for(i=0; i<nnodes; i++) {
		    fprintf(file_id, " %.12f %.12f %.12f\n", COORDS[i][0], COORDS[i][1], 0.0);
		}
		
		//Write coonectivity
		nsize = 5*nelem;
		fprintf(file_id, "CELLS %d %d\n", nelem, nsize);
		for(ie=0; ie<nelem; ie++) {
		    if( ELTYPE[ie] == VISIT_QUAD )
			fprintf(file_id, " %d %d %d %d %d\n", ie, MAPL2G[ie][0]-1, MAPL2G[ie][1]-1, MAPL2G[ie][2]-1, MAPL2G[ie][3]-1);
		    else if( ELTYPE[ie] == VISIT_TRIANGLE )
			fprintf(file_id, " %d %d %d %d\n", ie, MAPL2G[ie][0]-1, MAPL2G[ie][1]-1, MAPL2G[ie][2]-1);
		}
        
		fprintf(file_id, "CELL_TYPES %d\n", nelem);
		for(ie=0; ie<nelem; ie++)
		    fprintf(file_id, " %d\n", ELTYPE[ie]);
        
	    } else if( nsd == 3 ) {
		fprintf(file_id, "POINTS %d double\n", nnodes);
		for (i=0; i<nnodes; i++){
		    fprintf(file_id, " %.12f %.12f %.12f\n", COORDS[i][0], COORDS[i][1], COORDS[i][2]);
		}
            
		/*//COORDS --> COORDS1d:
		  j=1;
		  for (i=0; i<nnodes; i++){
                    
		  COORDS1d[j]   = COORDS[i][2];
		  COORDS1d[j+1] = COORDS[i][3];
		  COORDS1d[j+2] = COORDS[i][4];
		  j=j+3;
		  }
              
		  // Write the Coordinates 1D array to file
		  VECT2F_dVECT("COORDS1d.dat", COORDS1d, 1,nnodes*3);
		*/

		//Write coonectivity
		nsize = 9*nelem; //this is valid for hexa only
		fprintf(file_id, "CELLS %d %d\n", nelem, nsize);
		for(ie=0; ie<nelem; ie++) {
		    fprintf(file_id, " %d %d %d %d %d %d %d %d %d\n", 8, MAPL2G[ie][0]-1, MAPL2G[ie][1]-1, MAPL2G[ie][2]-1, MAPL2G[ie][3]-1, MAPL2G[ie][4]-1, MAPL2G[ie][5]-1, MAPL2G[ie][6]-1, MAPL2G[ie][7]-1);
		}
		
		fprintf(file_id, "CELL_TYPES %d\n", nelem);
		for(ie=0; ie<nelem; ie++)
		    fprintf(file_id, " %d\n", VISIT_HEXAHEDRON);
                
	    }
	} else {
	    /*
	     * Hight order grid by means of linear sub-elements
	     */
	    ncells = nelem*nop*nop*nop; //Total number of cells
	    fprintf(file_id, "POINTS %d double\n", nnodes);
	    for (i=0; i<nnodes; i++){
		fprintf(file_id, " %.12f %.12f %.12f\n", COORDS_HO[i][0], COORDS_HO[i][1], COORDS_HO[i][2]);
	    }
	   
	    nsize = 9*ncells;
	    fprintf(file_id, "CELLS %d %d\n", ncells, nsize);
	    for (int iel=0; iel<nelem; iel++) {
		for (int i=0; i<nop; i++) {
		    for (int j=0; j<nop; j++) {
			for (int k=0; k<nop; k++) {
			    
			    ii = min(i+1,nop-1);
			    jj = min(j+1,nop-1);
			    kk = min(k+1,nop-1);

			    fprintf(file_id, " 8 %d %d %d %d %d %d %d %d\n", \
				    CONN_HO2D[iel][ i][ j][ k],		\
				    CONN_HO2D[iel][ii][ j][ k],		\
				    CONN_HO2D[iel][ii][jj][ k],		\
				    CONN_HO2D[iel][ i][jj][ k],		\
				    CONN_HO2D[iel][ i][ j][kk],		\
				    CONN_HO2D[iel][ii][ j][kk],		\
				    CONN_HO2D[iel][ii][jj][kk],		\
				    CONN_HO2D[iel][ i][jj][kk]);
			    
			}
		    }
		}
	    }
	    
	}
    }
    fclose(file_id);
  
    return;
}

void dwrt2GMSH(char *file_name)
{  
    int i,j,ie,inode;
    int u,v,w;
    int elu,elv,elw, iel;
    int nsize, nfield_data;
  
    FILE *file_id;
  
    file_id = fopen(file_name, "w");
    if(file_id == NULL)
	printf("The file %s could not be open\n", file_name);
    else{
    
	//Open GMSH/ABAQUS file
	fprintf(file_id, "*Heading\n");
	fprintf(file_id, " %s\n", file_name);
    
	//Write coordinates
	fprintf(file_id, "*Node\n");
    
	//Write coordinates
	if( nsd == 2 )
	    {
		printf("\n ERROR (PRINT.c): GMSH now can only be written for nsd=3\n\n");
	    }
	else if( nsd == 3 )
	    {
		//fprintf(file_id, "$Nodes\n %d\n", nnodes);
		i=1;
	
		for (w=1; w<=nnodesz; w++){
		    for (v=1; v<=nnodesy; v++){
			for (u=1; u<=nnodesx; u++){
	    
			    fprintf(file_id, "%d, %.8f, %.8f, %.8f\n", i, COORDS[i][2], COORDS[i][3], COORDS[i][4]);
	      
			    i++;
			}
		    }
		}
	
		/*	for(i=1; i<=nnodes; i++)
			{
			fprintf(file_id, "%d %.12f %.12f %.12f\n",i, COORDS[i][2], COORDS[i][3], COORDS[i][4]);
			}
		*/
	
		//fprintf(file_id, "$EndNodes\n");
		//Write coonectivity
		if( ELTYPE[1] == VISIT_HEXAHEDRON )
		    nsize = 9*nelem;
		else if( ELTYPE[1] == VISIT_WEDGE )
		    nsize = 7*nelem;


		fprintf(file_id, "*Element, type=C3D8, ELSET=Volume1\n");
		for(ie=1; ie<=nelem; ie++)
		    {
			if( ELTYPE[ie-1] == VISIT_HEXAHEDRON )
			    fprintf(file_id, "%d, %d, %d, %d, %d, %d, %d, %d, %d\n", ie, MAPL2G[ie][2], MAPL2G[ie][3], MAPL2G[ie][4], MAPL2G[ie][5], MAPL2G[ie][6], MAPL2G[ie][7], MAPL2G[ie][8], MAPL2G[ie][9]);
			else if( ELTYPE[ie-1] == VISIT_WEDGE )
			    fprintf(file_id, " %d %d %d %d %d %d %d\n", 6, MAPL2G[ie][2], MAPL2G[ie][3], MAPL2G[ie][4], MAPL2G[ie][5], MAPL2G[ie][6], MAPL2G[ie][7]);
		    }
	
	    }
    }
    fclose(file_id);
  
    return;
}

void dwrt2CONN(char *file_name)
{  
    int i,j,ie,inode;
    int u,v,w;
    int elu,elv,elw, iel;
    int nsize, nfield_data;
  
    FILE *file_id;
  
    file_id = fopen(file_name, "w");
    if(file_id == NULL)
	printf("The file %s could not be open\n", file_name);
    else{
    
	fprintf(file_id, "%d\n", nelem);
    
	//Write coordinates
	if( nsd == 2 )
	    {
		printf("\n ERROR (PRINT.c): CONN now can only be written for nsd=3\n\n");
	    }
	else if( nsd == 3 )
	    {
		i=1;
	
		//Write coonectivity
		if( ELTYPE[1] == VISIT_HEXAHEDRON )
		    nsize = 9*nelem;
		else if( ELTYPE[1] == VISIT_WEDGE )
		    nsize = 7*nelem;
	
		for(ie=1; ie<=nelem; ie++)
		    {
			if( ELTYPE[ie-1] == VISIT_HEXAHEDRON )
			    fprintf(file_id, "%d %d %d %d %d %d %d %d\n", MAPL2G[ie][2], MAPL2G[ie][3], MAPL2G[ie][4], MAPL2G[ie][5], MAPL2G[ie][6], MAPL2G[ie][7], MAPL2G[ie][8], MAPL2G[ie][9]);
			else if( ELTYPE[ie-1] == VISIT_WEDGE )
			    fprintf(file_id, " %d %d %d %d %d %d\n", MAPL2G[ie][2], MAPL2G[ie][3], MAPL2G[ie][4], MAPL2G[ie][5], MAPL2G[ie][6], MAPL2G[ie][7]);
		    }
	
	    }
    }
    fclose(file_id);
  
    return;
}

int PRINT_WELCOME_MESSAGE(void)
{
    printf("\n #------------------------------------------------------------------#\n");
    printf(" # \n" );
    printf(" # Welcome to MMesh3D_M - V 2.0\n" );
    printf(" # \n" );
    printf(" # Author: Simone Marras, simone.marras@gmail.com\n" );
    printf(" #------------------------------------------------------------------#\n");

    return 0;
}

int PRINT_INFO(void)
{ 
    printf(" #------------------------------------------------------------------#\n #\n");
    printf(" # Input Entries from input file: %s\n", inputfile);
    if (lread_external_grid !=  0) {
	printf(" # Read GMSH external grid:   %s\n", external_grid_file_name);
    } else {
	printf(" # problem name:     %s\n", problem[0]);
	printf(" # meshing scheme:   %s %lf %lf\n", problem[1], parameters[7], parameters[8]);
	printf(" # element types:    %s\n", problem[2]);
	printf(" # nelx:             %d  %d\n", INPUTVariables[1], INPUTVariables[2]);
	printf(" # nely:             %d  %d\n", INPUTVariables[3], INPUTVariables[4]);
	printf(" # nelz:             %d  %d\n", INPUTVariables[5], INPUTVariables[6]);
	printf(" # nop:              %d\n", INPUTVariables[8]);
	printf(" # Periodicity x:    %s\n", problem[3]);
	printf(" # Periodicity y:    %s\n", problem[4]);
	printf(" # Periodicity z:    %s\n", problem[8]);
	printf(" # BDY NODES: (x,y,z)\n");
	for(i=1; i<=NBDY_NODES; i++)
	    printf(" #               %lf %lf %lf\n", BDY_COORDS[i][1],BDY_COORDS[i][2],BDY_COORDS[i][3]);
	printf(" # Topography:       %s %s %s\n", problem[9],problem[5],problem[11]);
	if( !strncmp(problem[9],"FUNCTION", 4) || !strncmp(problem[9],"function", 4) || !strncmp(problem[9],"Function", 4) || \
	    !strncmp(problem[9],"user", 4)     || !strncmp(problem[9],"USER", 4)     || !strncmp(problem[9],"User", 4))
	    {
		printf(" #   -Hmount:        %f\n", parameters[4]);
		printf(" #   -1/2 width x:   %f\n", parameters[5]);
		printf(" #   -1/2 width y:   %f\n", parameters[9]);
		printf(" #   -Lambda:        %f\n", parameters[6]);
	    }
	printf(" #\n");
    }
    printf(" # VTK OUTPUT:   \t %s\n", problem[7]);
    //printf(" # GMSH OUTPUT:  %s\n", problem[13]);
    //printf(" # BDY FILE:     %s\n", problem[12]);
    //printf(" # ALYA OUTPUT:  %s\n", problem[6]);
    printf(" #\n # End review of input entries \n");

    return 0;
}

void PRINT_ERROR(char *message)
{
    printf(" %s\n", message);
    exit(1);
}
