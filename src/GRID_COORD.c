/*GRID_COORD.c*/

/*Function that computes the grid coordinates

  simone.marras@gmail.com
  April 29, 2008
*/
#include <stdlib.h>
#include <string.h>

#include "GRID_COORD.h"
#include "GLOBAL_VARS.h"
#include "MYDEFINE.h"
#include "NRUTIL.h"


//void GRID_COORD(double **COORDS, int *grd_type, float xlength, float ylength, float zlength, float xmin_obs, float xmax_obs, float ymin_obs, float ymax_obs, float zmin_obs, float zmax_obs, int nnodesx, int nnodesy, int nnodesz, int iup, int nnodesx_obs, int idwn, int jup, int nnodesy_obs, int jdwn, int kup, int nnodesz_obs, int kdwn)
void GRID_COORD(void)
{
  //the INT GRD_TYPE of the function input defines the type of grid point distribution
  int i,j,k,count,ipoin;
  float deps;
	
  double *x, *y, *z;
  float xlength,ylength,zlength;

  //Bdy flag file: instead of sotirng the boundaries in an array, I write the flag 
  //to the file at every loop
  int BDYFLAG[4];
  FILE *BDYFLAGf_ID; 
  char *outfile_bdyflag;
  char nodesu[4], nodesv[4], nodesw[4], elorder[3]; //NOTE: if you want to store 3 digits, you need 4 spaces because C stores "/0" as last. If you don't, you will get the Abort Trap error!
  outfile_bdyflag = (char*) malloc(96 * sizeof(char *));
  
  //Bdy flag output file:
  sprintf(nodesu, "%d", nnodesx);
  sprintf(nodesv, "%d", nnodesy);
  sprintf(nodesw, "%d", nnodesz);
  sprintf(elorder,"%d", nop);
  
  strcpy(outfile_bdyflag, "3D_mesh_");
  strcat(outfile_bdyflag, problem[2]);
  strcat(outfile_bdyflag,"_");
  strcat(outfile_bdyflag,nodesu);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesv);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesw);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,elorder);
  strcat(outfile_bdyflag, "_BDYFLAG.dat");

  //Allocate:
  x = dvector(1,nnodesx);
  y = dvector(1,nnodesy);
  z = dvector(1,nnodesz);
  
  //WARNING: This will work correctly only for REGULAR, RECTANGULAR CUBES:
  xmin = BDY_COORDS[1][1];
  xmax = BDY_COORDS[2][1];
  xlength = xmax - xmin;
  ymin = BDY_COORDS[1][2];
  ymax = BDY_COORDS[4][2];
  ylength = ymax - ymin;
  zmin = BDY_COORDS[1][3];
  zmax = BDY_COORDS[5][3];
  zlength = zmax - zmin;

  /*COMPUTATIONS*/
  //A: NO OBSTACLE in the Flow:
  x[1] = 0;
  dx = xlength/(nnodesx-1);
  //Grid coordinates
  for(i=2; i<=nnodesx; i++){
    dx = xlength/(nnodesx-1);
    x[i] = x[i-1] + dx;
    //printf("dx[%d] = %lf\t x[%d] = %lf\n",i, dx[i], i, x[i]);
  }
  y[1] = 0;
  //printf(" ymin = %f\n ymax = %f\n", ymin, ymax);
  dy = ylength/(nnodesy-1);
  //Grid coordinates
  for(i=2; i<=nnodesy; i++){
    dy = ylength/(nnodesy-1);
    y[i] = y[i-1] + dy;
  }
  z[1] = 0;
  dz = zlength/(nnodesz-1);
  //Grid coordinates
  for(i=2; i<=nnodesz; i++){
    dz = zlength/(nnodesz-1);
    z[i] = z[i-1] + dz;
  }
  
  /*****************************************************************************
   * Open the bdy flag file for writing:
   *****************************************************************************/
  BDYFLAGf_ID = fopen(outfile_bdyflag, "a+");
  
  /*****************************************************************************
   * Mesh the domain:
   *****************************************************************************/
  
  //Build the COORDS[1:npoin][1:nsd] array
  ipoin = 0;
  for(k=1; k<=nnodesz; k++){
    for(j=1; j<=nnodesy; j++){
      for(i=1; i<=nnodesx; i++){     
	
	ipoin = ipoin + 1;
	
	COORDS[ipoin][1] = ipoin;
	COORDS[ipoin][2] = x[i];
	COORDS[ipoin][3] = y[j];
	COORDS[ipoin][4] = z[k];

	//Assign the boundary flags:
	//a) BOTTOM CORNERS:
	if( i == 1 && j == 1 && k == 1){
	  BDYFLAG[0] = 1;				// == 1 means that is it a boundary node on the front-bott-left corner.
	}
	else if( i == nnodesx && j == 1 && k == 1){
	  BDYFLAG[0] = 2;				// == 2 means that is it a boundary node on the front-bott-right corner.
	}
	else if( i == nnodesx && j == nnodesy && k == 1){
	  BDYFLAG[0] = 3;				// == 3 means that is it a boundary node on the back-bott-right corner.
	}
	else if( i == 1 && j == nnodesy && k == 1){
	  BDYFLAG[0] = 4;				// == 4 means that is it a boundary node on the back-bott-left corner.
	}
	//b) TOP CORNERS:
	else if( i == 1 && j == 1 && k == nnodesz){
	  BDYFLAG[0] = 8;				// == 8 means that is it a boundary node on the front-top-left corner.
	}
	else if( i == nnodesx && j == 1 && k == nnodesz){
	  BDYFLAG[0] = 5;				// == 5 means that is it a boundary node on the front-top-right corner.
	}
	else if( i == nnodesx && j == nnodesy && k == nnodesz){
	  BDYFLAG[0] = 6;				// == 6 means that is it a boundary node on the back-top-right corner.
	}
	else if( i == 1 && j == nnodesy && k == nnodesz){
	  BDYFLAG[0] = 7;				// == 7 means that is it a boundary node on the back-top-left corner.
	}
	
	//c) BOTT-EDGES (horizontal):
	else if( j == 1 && k == 1 && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 9;				// == 9 means that is it a boundary node on the bott-front edge.
	}
	else if( i == nnodesx && k == 1 && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 10;				// == 10 means that is it a boundary node on the bottom-right edge.
	}
	else if( j == nnodesy && k == 1 && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 11;				// == 11 means that is it a boundary node on the bottom-back edge.
	}
	else if( i == 1 && k == 1 && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 12;				// == 12 means that is it a boundary node on the bottom-left edge.
	}
	
	//d) TOP-EDGES (horizontal):
	else if( j == 1 && k == nnodesz && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 17;				// == 17 means that is it a boundary node on the top-front edge.
	}
	else if( i == nnodesx && k == nnodesz && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 18;				// == 18 means that is it a boundary node on the top-right edge.
	}
	else if( j == nnodesy && k == nnodesz && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 19;				// == 19 means that is it a boundary node on the top-back edge.
	}
	else if( i == 1 && k == nnodesz && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 20;				// == 20 means that is it a boundary node on the top-left edge.
	}
	
	//e) JERTICAL-EDGES:
	else if( i == 1 && j == 1 && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 16;				// == 16 means that is it a boundary node on the front-left jedge.
	}
	else if( i == nnodesx && j == 1 && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 13;				// == 13 means that is it a boundary node on the front-right jedge.
	}
	else if( i == nnodesx && j == nnodesy && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 14;				// == 14 means that is it a boundary node on the back-right jedge.
	}
	else if( i == 1 && j == nnodesy && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 15;				// == 15 means that is it a boundary node on the back-left jedge.
	}
	
	//f) FACES:
	else if( k == 1 && (i<nnodesx && i>1) && (j<nnodesy && j>1)){
	  BDYFLAG[0] = 21;				// == 21 means that is it a boundary node on the bottom face
	}
	else if( k == nnodesz && (i<nnodesx && i>1) && (j<nnodesy && j>1)){
	  BDYFLAG[0] = 23;				// == 23 means that is it a boundary node on the top face
	}
	else if( i == nnodesx && (j>1 && j<nnodesy) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 22;				// == 22 means that is it a boundary node on the east face
	}
	else if( i == 1 && (j>1 && j<nnodesy) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 24;				// == 24 means that is it a boundary node on the kest face
	}
	else if( j == 1 && (i>1 && i<nnodesx) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 25;				// == 25 means that is it a boindary node on the front face
	}
	else if( j == nnodesy && (i>1 && i<nnodesx) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 26;				// == 26 means that is it a boindary node on the back face
	}
	else
	  {
	    BDYFLAG[0] = 0;		                // == 0 means that is it NOT a boundary node
	  }
	//End assign bdy flags.
	fprintf(BDYFLAGf_ID, "%d %lf %lf %lf\n", BDYFLAG[0], COORDS[ipoin][2],COORDS[ipoin][3],COORDS[ipoin][4]);
	
      }
    }
  }
  
  //Close the bdy file after use:
  fclose(BDYFLAGf_ID);
  
  /*****************************/
  /*END Coordinates computation*/
  /*****************************/
  
  //Free allocation:
  free_dvector(x, 1,nnodesx);
  free_dvector(y, 1,nnodesy);
  free_dvector(z, 1,nnodesz);
  free(outfile_bdyflag);

  return;
}


/******************************/
/*NON staggered grid coords****/
/******************************//*
 * HIGH-ORDER
 */
void GRID_COORD_HIGH_ORDER(void)
{
  //the INT GRD_TYPE of the function input defines the type of grid point distribution
  int i,j,k,count,ipoin;
  int ngl;
  float deps;
	
  double *x, *y, *z;
  float xlength,ylength,zlength;

  //Bdy flag file: instead of sotirng the boundaries in an array, I write the flag 
  //to the file at every loop
  int BDYFLAG[4];
  FILE *BDYFLAGf_ID; 
  char *outfile_bdyflag;
  char nodesu[4];
  char nodesv[4];
  char nodesw[4];
  char elorder[3]; //NOTE: if you want to store 3 digits, you need 4 spaces because C stores "/0" as last. If you don't, you will get the Abort Trap error!
  
  //Allocate string for bdy output file:
  outfile_bdyflag = (char*) malloc(96 * sizeof(char *));

  //Local
  int ix, iy, iz, ip, ib, ip1, ip2, ies;
  int ielex, ieley, ielez;
  int iglx, igly, iglz;
  int ilglz, ilgly, ilglx,ie;
  int i1, j1, k1, ii, jj, kk, ih;
  int nglm1;
  
  //Bdy flag output file:
  sprintf(nodesu, "%d", nnodesx);
  sprintf(nodesv, "%d", nnodesy);
  sprintf(nodesw, "%d", nnodesz);
  sprintf(elorder,"%d", nop);
  
  strcpy(outfile_bdyflag, "3D_mesh_");
  strcat(outfile_bdyflag, problem[2]);
  strcat(outfile_bdyflag,"_");
  strcat(outfile_bdyflag,nodesu);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesv);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesw);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,elorder);
  strcat(outfile_bdyflag, "_BDYFLAG.dat");

  //Allocate:
  x = dvector(1,nnodesx);
  y = dvector(1,nnodesy);
  z = dvector(1,nnodesz);
  
  //WARNING: This will work correctly only for REGULAR, RECTANGULAR CUBES:
  xmin = BDY_COORDS[1][1];
  xmax = BDY_COORDS[2][1];
  xlength = xmax - xmin;
  ymin = BDY_COORDS[1][2];
  ymax = BDY_COORDS[4][2];
  ylength = ymax - ymin;
  zmin = BDY_COORDS[1][3];
  zmax = BDY_COORDS[5][3];
  zlength = zmax - zmin;

  /*COMPUTATIONS*/
  //A: NO OBSTACLE in the Flow:
  x[1] = 0;
  dx = xlength/(nnodesx-1);
  //Grid coordinates
  for(i=2; i<=nnodesx; i++){
    dx = xlength/(nnodesx-1);
    x[i] = x[i-1] + dx;
    //printf("dx[%d] = %lf\t x[%d] = %lf\n",i, dx[i], i, x[i]);
  }
  y[1] = 0;
  //printf(" ymin = %f\n ymax = %f\n", ymin, ymax);
  dy = ylength/(nnodesy-1);
  //Grid coordinates
  for(i=2; i<=nnodesy; i++){
    dy = ylength/(nnodesy-1);
    y[i] = y[i-1] + dy;
  }
  z[1] = 0;
  dz = zlength/(nnodesz-1);
  //Grid coordinates
  for(i=2; i<=nnodesz; i++){
    dz = zlength/(nnodesz-1);
    z[i] = z[i-1] + dz;
  }
  
  /*****************************************************************************
   * Open the bdy flag file for writing:
   *****************************************************************************/
  BDYFLAGf_ID = fopen(outfile_bdyflag, "a+");
  
  /*****************************************************************************
   * Mesh the domain:
   *****************************************************************************/
  
  ngl   = nop + 1;
  nglm1 = ngl - 1;

  // Generate Coordinates Here 
  ip = 1; 
  ii = 0;
  jj = 0;
  kk = 0;
  
  //Build the COORDS[1:npoin][1:nsd] array
  ipoin = 0;
  for(k=1; k<=nnodesz; k++){
    for(j=1; j<=nnodesy; j++){
      for(i=1; i<=nnodesx; i++){     
	
	ipoin = ipoin + 1;
	
	COORDS[ipoin][1] = ipoin;
	COORDS[ipoin][2] = x[i];
	COORDS[ipoin][3] = y[j];
	COORDS[ipoin][4] = z[k];

	//Assign the boundary flags:
	//a) BOTTOM CORNERS:
	if( i == 1 && j == 1 && k == 1){
	  BDYFLAG[0] = 1;				// == 1 means that is it a boundary node on the front-bott-left corner.
	}
	else if( i == nnodesx && j == 1 && k == 1){
	  BDYFLAG[0] = 2;				// == 2 means that is it a boundary node on the front-bott-right corner.
	}
	else if( i == nnodesx && j == nnodesy && k == 1){
	  BDYFLAG[0] = 3;				// == 3 means that is it a boundary node on the back-bott-right corner.
	}
	else if( i == 1 && j == nnodesy && k == 1){
	  BDYFLAG[0] = 4;				// == 4 means that is it a boundary node on the back-bott-left corner.
	}
	//b) TOP CORNERS:
	else if( i == 1 && j == 1 && k == nnodesz){
	  BDYFLAG[0] = 8;				// == 8 means that is it a boundary node on the front-top-left corner.
	}
	else if( i == nnodesx && j == 1 && k == nnodesz){
	  BDYFLAG[0] = 5;				// == 5 means that is it a boundary node on the front-top-right corner.
	}
	else if( i == nnodesx && j == nnodesy && k == nnodesz){
	  BDYFLAG[0] = 6;				// == 6 means that is it a boundary node on the back-top-right corner.
	}
	else if( i == 1 && j == nnodesy && k == nnodesz){
	  BDYFLAG[0] = 7;				// == 7 means that is it a boundary node on the back-top-left corner.
	}
	
	//c) BOTT-EDGES (horizontal):
	else if( j == 1 && k == 1 && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 9;				// == 9 means that is it a boundary node on the bott-front edge.
	}
	else if( i == nnodesx && k == 1 && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 10;				// == 10 means that is it a boundary node on the bottom-right edge.
	}
	else if( j == nnodesy && k == 1 && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 11;				// == 11 means that is it a boundary node on the bottom-back edge.
	}
	else if( i == 1 && k == 1 && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 12;				// == 12 means that is it a boundary node on the bottom-left edge.
	}
	
	//d) TOP-EDGES (horizontal):
	else if( j == 1 && k == nnodesz && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 17;				// == 17 means that is it a boundary node on the top-front edge.
	}
	else if( i == nnodesx && k == nnodesz && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 18;				// == 18 means that is it a boundary node on the top-right edge.
	}
	else if( j == nnodesy && k == nnodesz && (i>1 && i<nnodesx)){
	  BDYFLAG[0] = 19;				// == 19 means that is it a boundary node on the top-back edge.
	}
	else if( i == 1 && k == nnodesz && (j>1 && j<nnodesy)){
	  BDYFLAG[0] = 20;				// == 20 means that is it a boundary node on the top-left edge.
	}
	
	//e) JERTICAL-EDGES:
	else if( i == 1 && j == 1 && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 16;				// == 16 means that is it a boundary node on the front-left jedge.
	}
	else if( i == nnodesx && j == 1 && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 13;				// == 13 means that is it a boundary node on the front-right jedge.
	}
	else if( i == nnodesx && j == nnodesy && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 14;				// == 14 means that is it a boundary node on the back-right jedge.
	}
	else if( i == 1 && j == nnodesy && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 15;				// == 15 means that is it a boundary node on the back-left jedge.
	}
	
	//f) FACES:
	else if( k == 1 && (i<nnodesx && i>1) && (j<nnodesy && j>1)){
	  BDYFLAG[0] = 21;				// == 21 means that is it a boundary node on the bottom face
	}
	else if( k == nnodesz && (i<nnodesx && i>1) && (j<nnodesy && j>1)){
	  BDYFLAG[0] = 23;				// == 23 means that is it a boundary node on the top face
	}
	else if( i == nnodesx && (j>1 && j<nnodesy) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 22;				// == 22 means that is it a boundary node on the east face
	}
	else if( i == 1 && (j>1 && j<nnodesy) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 24;				// == 24 means that is it a boundary node on the kest face
	}
	else if( j == 1 && (i>1 && i<nnodesx) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 25;				// == 25 means that is it a boindary node on the front face
	}
	else if( j == nnodesy && (i>1 && i<nnodesx) && (k>1 && k<nnodesz)){
	  BDYFLAG[0] = 26;				// == 26 means that is it a boindary node on the back face
	}
	else
	  {
	    BDYFLAG[0] = 0;		                // == 0 means that is it NOT a boundary node
	  }
	//End assign bdy flags.
	fprintf(BDYFLAGf_ID, "%d %lf %lf %lf\n", BDYFLAG[0], COORDS[ipoin][2],COORDS[ipoin][3],COORDS[ipoin][4]);
	
      }
    }
  }

  //Close the bdy file after use:
  fclose(BDYFLAGf_ID);
  
  /*****************************/
  /*END Coordinates computation*/
  /*****************************/
  
  //Free allocation:
  free_dvector(x, 1,nnodesx);
  free_dvector(y, 1,nnodesy);
  free_dvector(z, 1,nnodesz);
  free(outfile_bdyflag);

  return;
}/* END high-order */
