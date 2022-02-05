/***************************************
 *
 *	function CONN = grid2conn(nnodesx, nnodesy)
 *
 * This function builds the connectivity matrix for a 2D structured grid of QUADs or TRIs
 * simone.marras@gmail.com
 * 
 * 2009
 */
#include <stdio.h>
#include "GLOBAL_VARS.h"
#include "NRUTIL.h"

/* global variables */
int kstrt, jstrt, istrt;

#define TRI   3
#define QUAD  4
#define HEXA  8
#define WEDGE 6

#define INDX_OFFSET 1

/********************************************************************************************
 *
 * HEXAs:
 *
 ********************************************************************************************/
int GRID2CONNhex(unsigned int strting_elem_number, int *ELTYPE, int iblock)
{
  int i,j,k;
	
  int nelx = nnodesx-1;
  int nely = nnodesy-1;
  int nelz = nnodesz-1;
  int nelem = nelx*nely*nelz;
  int flag,counter,index_end,iel;

  nbdy_faces = 2*nelx*nelz + 2*nely*nelz + 2*nelx*nelz;
  
  if( iblock == 1){
    istrt = 0;
    jstrt = 0;
    kstrt = 0;
  }
  
  //Initialize:
  for (i=0; i<nelem; i++){
    ELTYPE[i] = HEXA;
    for (j=0; j<8; j++){
	CONN[i][j] = 0;
    }
  }
  //END Initialize
  flag = 0;
  if( iblock > 1)
    flag = 1;
  
  kstrt = 0;
  jstrt = 0;
  index_end = 0;
  iel = 0;
  
  //Loop to defined the connectivity
  //Base layer only (k = 1):
  k = 1;
  for (j=1+jstrt; j<=nnodesy-1+jstrt; j++){
    for (i=1; i<=nnodesx-1; i++){
      
      //CONN[iel][0] = iel;
      //CONN[iel][0] = iel + (j - 1)+1;
      
      CONN[iel][0] =  iel + index_end + (j - 1)+1;
      CONN[iel][1] = CONN[iel][0] + 1;
      CONN[iel][2] = CONN[iel][1] + nnodesx;
      CONN[iel][3] = CONN[iel][2] - 1;
      
      CONN[iel][4] = CONN[iel][0] + nnodesx*nnodesy;
      CONN[iel][5] = CONN[iel][1] + nnodesx*nnodesy;
      CONN[iel][6] = CONN[iel][2] + nnodesx*nnodesy;
      CONN[iel][7] = CONN[iel][3] + nnodesx*nnodesy;

      iel = iel + 1;
    }
  }
  strting_elem_number = iel;
  index_end = nnodesx*nnodesy + 1;

  //Loop to define the connectivity
  //Layer k > 1:
  for (k=2; k<=nnodesz-1; k++){
    for (j=1; j<=nnodesy-1; j++){
      for (i=1; i<=nnodesx-1; i++){
	
	//CONN[iel][0] = iel;
	
	CONN[iel][0] = CONN[iel - nelx*nely][4];   //It uses the info of the top face of the undelrying element.
	CONN[iel][1] = CONN[iel][0] + 1;
	CONN[iel][2] = CONN[iel][1] + nnodesx;
	CONN[iel][3] = CONN[iel][2] - 1;
	
	CONN[iel][4] = CONN[iel][0] + nnodesx*nnodesy;
	CONN[iel][5] = CONN[iel][1] + nnodesx*nnodesy;
	CONN[iel][6] = CONN[iel][2] + nnodesx*nnodesy;
	CONN[iel][7] = CONN[iel][3] + nnodesx*nnodesy;

	iel = iel + 1;
	index_end = index_end + 1;
      }
    }
    index_end = index_end + nnodesx*nnodesy;
  }

  if(nop < 2) {
      for (int iel=0; iel<nelem; iel++) {
	  for (int i=0; i<HEXA; i++) {
	      MAPL2G[iel][i] = CONN[iel][i];
	  }
      }
  }
  
  if( iblock == 1){
    jstrt = j-1; //Stored from the previous block
    kstrt = k-1;
  }
  
  return 0;
}

