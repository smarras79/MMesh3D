/***************************************
 *
 *	function CONN = grid2conn(nnodesx, nnodesy)
 *
 * This function builds the connectivity matrix for a 2D structured grid of QUADs or TRIs
 * simone.marras@gmail.com
 */

#include <stdio.h>

#include "NRUTIL.h"

/* global variables */
int kstrt, jstrt, istrt;

#define TRI   3
#define QUAD  4
#define HEXA  8
#define WEDGE 6

/********************************************************************************************
 *
 * QUADS
 *
 ********************************************************************************************/
int GRID2CONN(unsigned int strting_node_number, unsigned int strting_elem_number, 
	      int nnodesx, int nnodesy, int **CONN, int *ELTYPE, int iblock)
{
  //int **CONN;
	
  int i,j,k;
	
  int nelx = nnodesx-1;
  int nely = nnodesy-1;
  int nelem = nelx*nely;
  int flag;
	
  //printf(" nnodes = %d nelem = %d strting_elem_number = %d\n", nnodesx*nnodesy, nelem, strting_elem_number);

	
  if( iblock == 1){
    jstrt = 0;
  }
	
  //Initialize:
  for (i = strting_elem_number+1; i<=nelem+strting_elem_number; i++){
    for (j = 1; j<=5; j++){
      CONN[i][j] = 0;
    }
  }
  //END Initialize	

  flag = 0;
  if( iblock > 1)
    flag = 1;
	
  int iel = strting_elem_number+1;
  for (j = 1+jstrt; j<=nnodesy-1+jstrt; j++){
    for (i = 1; i<=nnodesx-1; i++){

      //printf(" j = %d  i = %d iel = %d\n", j, i, iel);
			
      CONN[iel][1] = iel;
      CONN[iel][2] = iel + (j - 1);
      CONN[iel][3] = CONN[iel][2] + 1;
      CONN[iel][4] = CONN[iel][3] + nnodesx;
      CONN[iel][5] = CONN[iel][4] - 1;
			
      //Modify this if you add the option of creating triangular elements:
      ELTYPE[iel] = 4;
			
      iel = iel + 1;
	
    }
  }

  if( iblock == 1){
    jstrt = j-1; //Stored from the previous block
  }

  return 0;
}


/********************************************************************************************
 *
 * TRIs
 * Numbering algorithm for triangular elements obtained from structured quads:
 *
 *	[1] Date A., "Introduction to Computational Fluid Dynamics", Cambridge Uni. Press. 1st Ed.
 *				 Ch. 8.5, p. 250.
 ********************************************************************************************/
int GRID2CONNtri(unsigned int strting_node_number, unsigned int strting_elem_number, 
		 int nnodesx, int nnodesy, int nelem, int **CONN, int *ELTYPE, int iblock)
{
  //int **CONN;
	
  int i,j,k;
  unsigned int problem_type = 0;
	
  int NE1, NE2, NV, NV1, NV2, NV3, NV4;
  int M;
	
  if( iblock == 1){
    jstrt = 0;
  }
	
  //Ask the user whether the TRIs are all with the same
  //inclination or if they are numbered differently as a function
  //of the node of reference:
  //
  //	if the counting node of element i-th is ODD then the original QUAD
  //	is divided with a diagonal going from bottom-left to top-right.
  //	Otherwise, if the node of elem. i-th is EVEN then the original QUAD
  //	is divided with a diagonal going from bottom-right to top-left.
  //
  //	See [1] for details.
  //
	
  printf("\n # USER INPUTS on the TRI types: \n");
  printf(" # Select 0 or 1:\n ");
  printf(" # 0 -> Default (The diagonal is decided by the ODD or EVEN ref. node)\n ");
  printf(" # 1 -> All elements are of the same geometry\n ");
  //printf(" # 2 -> The diagonal is decided by the ODD or EVEN character of the node of reference\n\n ");
  scanf("%d", &problem_type); //NOTA BENE: "%d" NON DEVE AVER NESSUNO SPAZIO TRA LA 'd' e le " di chiusura o il codice va in palla!


	
  //Initialize:
  for (i = strting_elem_number+1; i<=nelem+strting_elem_number; i++){
    for (j = 1; j<=4; j++){
      CONN[i][j] = 0;
    }
  }
  //END Initialize	
	
  //NE1 = 0;
  NE1 = strting_elem_number;
  for (j = 1+jstrt; j<=nnodesy-1+jstrt; j++){
    for (i = 1; i<=nnodesx-1; i++){
			
      NV = i + (j-1)*nnodesx;
      NE1 = NE1 + 1;
      NE2 = NE1 + 1;
			
      //printf(" NE1 =%d,  NE2=%d\n", NE1, NE2);
			
      if(problem_type == 1)
	M = 3;
      else if(problem_type == 0)    //(now the user input 0 or 2 gives the same selection of M)
	M = NV - ((int)(NV/2))*2; //Remainder of NV/2
			
      NV1 = NV;
      NV2 = NV1 + nnodesx;
      NV3 = NV2 + 1;
      NV4 = NV1 + 1;
			
      if( M == 1) //Meaning that the node number NV is odd
	{
	  CONN[NE1][1] = NE1;
	  CONN[NE1][2] = NV1;
	  CONN[NE1][3] = NV3;
	  CONN[NE1][4] = NV2;
				
	  CONN[NE2][1] = NE2;
	  CONN[NE2][2] = NV1;
	  CONN[NE2][3] = NV4;
	  CONN[NE2][4] = NV3;
				
	}
      else if(M == 0) //Meaning that the node number NV is even
	{
	  CONN[NE1][1] = NE1;
	  CONN[NE1][2] = NV1;
	  CONN[NE1][3] = NV4;
	  CONN[NE1][4] = NV2;
				
	  CONN[NE2][1] = NE2;
	  CONN[NE2][2] = NV4;
	  CONN[NE2][3] = NV3;
	  CONN[NE2][4] = NV2;
			
	}
      else{
	CONN[NE1][1] = NE1;
	CONN[NE1][2] = NV1;
	CONN[NE1][3] = NV3;
	CONN[NE1][4] = NV2;
				
	CONN[NE2][1] = NE2;
	CONN[NE2][2] = NV1;
	CONN[NE2][3] = NV4;
	CONN[NE2][4] = NV3;
      }
			
      //Modify this if you add the option of creating triangular elements:
      ELTYPE[NE1] = 3;
      ELTYPE[NE2] = 3;
							
      NE1 = NE2;

    }
  }
	
  if( iblock == 1){
    jstrt = j-1; //Stored from the previous block
  }

  return 0;
}


/********************************************************************************************
 *
 * HEXAHEDRONS:
 *
 ********************************************************************************************/
int GRID2CONNhex(unsigned int strting_node_number, unsigned int strting_elem_number, 
		 int nnodesx, int nnodesy, int nnodesz, int **CONN, int *ELTYPE, int iblock)
{
  int i,j,k;
	
  int nelx = nnodesx-1;
  int nely = nnodesy-1;
  int nelz = nnodesz-1;
  int nelem = nelx*nely*nelz;
  int flag,counter,index_end,iel;
	
  //printf(" nnodes = %d nelem = %d strting_elem_number = %d\n", nnodesx*nnodesy, nelem, strting_elem_number);
  //printf(" GRID2CONNhex: nx %d ny %d\n", nnodesx, nnodesy);

  if( iblock == 1){
    istrt = 0;
    jstrt = 0;
    kstrt = 0;
  }

  //Initialize:
  for (i = strting_elem_number+1; i<=nelem+strting_elem_number; i++){
    ELTYPE[i] = HEXA;
    for (j = 1; j<=HEXA+1; j++){
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
  iel = strting_elem_number;
  
  //Loop to defined the connectivity
  //Base layer only (k = 1):
  k = 1;
  for (j=1+jstrt; j<=nnodesy-1+jstrt; j++){
    for (i=1; i<=nnodesx-1; i++){
      
      iel = iel + 1;
      CONN[iel][1] = iel;
      CONN[iel][2] = iel + (j - 1);
      
      CONN[iel][2] =  iel + index_end + (j - 1);
      CONN[iel][3] = CONN[iel][2] + 1;
      CONN[iel][4] = CONN[iel][3] + nnodesx;
      CONN[iel][5] = CONN[iel][4] - 1;
      
      CONN[iel][6] = CONN[iel][2] + nnodesx*nnodesy;
      CONN[iel][7] = CONN[iel][3] + nnodesx*nnodesy;
      CONN[iel][8] = CONN[iel][4] + nnodesx*nnodesy;
      CONN[iel][9] = CONN[iel][5] + nnodesx*nnodesy;
    }
  }
  strting_elem_number = iel;
  index_end = nnodesx*nnodesy + 1;

  //Loop to define the connectivity
  //Layer k > 1:
  for (k=2; k<=nnodesz-1; k++){
    for (j=1; j<=nnodesy-1; j++){
      for (i=1; i<=nnodesx-1; i++){
	
	iel = iel + 1;
	CONN[iel][1] = iel;
	
	CONN[iel][2] = CONN[iel - nelx*nely][6];   //It uses the info of the top face of the undelrying element.
	CONN[iel][3] = CONN[iel][2] + 1;
	CONN[iel][4] = CONN[iel][3] + nnodesx;
	CONN[iel][5] = CONN[iel][4] - 1;
	
	CONN[iel][6] = CONN[iel][2] + nnodesx*nnodesy;
	CONN[iel][7] = CONN[iel][3] + nnodesx*nnodesy;
	CONN[iel][8] = CONN[iel][4] + nnodesx*nnodesy;
	CONN[iel][9] = CONN[iel][5] + nnodesx*nnodesy;

	index_end = index_end + 1;
      }
    }
    index_end = index_end + nnodesx*nnodesy;
  }
  
  if( iblock == 1){
    jstrt = j-1; //Stored from the previous block
    kstrt = k-1;
  }
  
  return 0;
}


/********************************************************************************************
 *
 * WEDGEs with triangular base
 *
 ********************************************************************************************/
int GRID2CONNwedge(unsigned int strting_node_number, unsigned int strting_elem_number, 
		   int nnodesx, int nnodesy, int nnodesz, int **CONN, int *ELTYPE, int iblock)
{
  int i,j,k,iel;
  int flag,counter,index_end;
  unsigned int problem_type = 0;
  
  int nelx = nnodesx-1;
  int nely = nnodesy-1;
  int nelz = nnodesz-1;
  int nelem = 2*nelx*nely*nelz;

  int NE1, NE2, NV, NV1, NV2, NV3, NV4;
  int M;
	
  if( iblock == 1){
    jstrt = 0;
  }
	
  //Ask the user whether the TRIs are all with the same
  //inclination or if they are numbered differently as a function
  //of the node of reference:
  //
  //	if the counting node of element i-th is ODD then the original QUAD
  //	is divided with a diagonal going from bottom-left to top-right.
  //	Otherwise, if the node of elem. i-th is EVEN then the original QUAD
  //	is divided with a diagonal going from bottom-right to top-left.
  //
  //	See [1] for details.
  //

  printf("\n # USER INPUTS on the TRI types: \n");
  printf(" # Select 0 or 1:\n ");
  printf(" # 0 -> Default (The diagonal is decided by the ODD or EVEN ref. node)\n ");
  printf(" # 1 -> All elements are of the same geometry\n ");
  //printf(" # 2 -> The diagonal is decided by the ODD or EVEN character of the node of reference\n\n ");
  scanf("%d", &problem_type); 
  //NOTA BENE: "%d" NON DEVE AVER NESSUNO SPAZIO TRA LA 'd' e le " di chiusura o il codice va in palla!
 
  //Initialize:
  for (i = strting_elem_number+1; i<=nelem+strting_elem_number; i++){
    //    ELTYPE[i] = WEDGE;
    for (j = 1; j<=WEDGE+1; j++){
      CONN[i][j] = 1;
    }
  }
  //END Initialize
  
  //Base layer only (k = 1):
  //NE1 = 0;
  iel = strting_elem_number;
  NE1 = strting_elem_number;
  for (j = 1+jstrt; j<=nnodesy-1+jstrt; j++){
    for (i = 1; i<=nnodesx-1; i++){
			
      NV = i + (j-1)*nnodesx;
      NE1 = NE1 + 1;
      NE2 = NE1 + 1;
			
      iel = iel + 1;
			
      if(problem_type == 1)
	M = 3;
      else if(problem_type == 0)    //(now the user input 0 or 2 gives the same selection of M)
	M = NV - ((int)(NV/2))*2; //Remainder of NV/2
			
      NV1 = NV;
      NV2 = NV1 + nnodesx;
      NV3 = NV2 + 1;
      NV4 = NV1 + 1;
      
      if( M == 1) //Meaning that the node number NV is odd
	{
	  CONN[NE1][1] = NE1;
	
	  CONN[NE1][2] = NV1;
	  CONN[NE1][3] = NV3;
	  CONN[NE1][4] = NV2;
	
	  CONN[NE1][5] = CONN[NE1][2] + nnodesx*nnodesy;
	  CONN[NE1][6] = CONN[NE1][3] + nnodesx*nnodesy;
	  CONN[NE1][7] = CONN[NE1][4] + nnodesx*nnodesy;
	  	  
	  CONN[NE2][1] = NE2;
	 
	  CONN[NE2][2] = NV1;
	  CONN[NE2][3] = NV4;
	  CONN[NE2][4] = NV3;

	  CONN[NE2][5] = CONN[NE2][2] + nnodesx*nnodesy;
	  CONN[NE2][6] = CONN[NE2][3] + nnodesx*nnodesy;
	  CONN[NE2][7] = CONN[NE2][4] + nnodesx*nnodesy;
	 	  
	}
      else if(M == 0) //Meaning that the node number NV is even
	{
	  CONN[NE1][1] = NE1;

	  CONN[NE1][2] = NV1;
	  CONN[NE1][3] = NV4;
	  CONN[NE1][4] = NV2;

	  CONN[NE1][5] = CONN[NE1][2] + nnodesx*nnodesy;
	  CONN[NE1][6] = CONN[NE1][3] + nnodesx*nnodesy;
	  CONN[NE1][7] = CONN[NE1][4] + nnodesx*nnodesy;
	  
	  CONN[NE2][1] = NE2;

	  CONN[NE2][2] = NV4;
	  CONN[NE2][3] = NV3;
	  CONN[NE2][4] = NV2;

	  CONN[NE2][5] = CONN[NE2][2] + nnodesx*nnodesy;
	  CONN[NE2][6] = CONN[NE2][3] + nnodesx*nnodesy;
	  CONN[NE2][7] = CONN[NE2][4] + nnodesx*nnodesy;
	  
	}
      else{
	
	CONN[NE1][1] = NE1;

	CONN[NE1][2] = NV1;
	CONN[NE1][3] = NV3;
	CONN[NE1][4] = NV2;

	CONN[NE1][5] = CONN[NE1][2] + nnodesx*nnodesy;
	CONN[NE1][6] = CONN[NE1][3] + nnodesx*nnodesy;
	CONN[NE1][7] = CONN[NE1][4] + nnodesx*nnodesy;

	CONN[NE2][1] = NE2;

	CONN[NE2][2] = NV1;
	CONN[NE2][3] = NV4;
	CONN[NE2][4] = NV3;
	
	CONN[NE2][5] = CONN[NE2][2] + nnodesx*nnodesy;
	CONN[NE2][6] = CONN[NE2][3] + nnodesx*nnodesy;
	CONN[NE2][7] = CONN[NE2][4] + nnodesx*nnodesy;
      }
			
      //Modify this if you add the option of creating triangular elements:
      ELTYPE[NE1] = WEDGE;
      ELTYPE[NE2] = WEDGE;
							
      NE1 = NE2;

    }
  }

  strting_elem_number = iel;
  index_end = nnodesx*nnodesy + 1;
 
    //Loop to define the connectivity
  //Layer k > 1:
  for (k=2; k<=nnodesz-1; k++){
    for (j=1; j<=nnodesy-1; j++){
      for (i=1; i<=nnodesx-1; i++){

	NE1 = NE1 + 1;
	NE2 = NE1 + 1;

	CONN[NE1][1] = NE1;
	CONN[NE1][2] = CONN[NE1 - 2*nelx*nely][5];
	CONN[NE1][3] = CONN[NE1 - 2*nelx*nely][6];
	CONN[NE1][4] = CONN[NE1 - 2*nelx*nely][7];

	CONN[NE1][5] = CONN[NE1][2] + nnodesx*nnodesy;
	CONN[NE1][6] = CONN[NE1][3] + nnodesx*nnodesy;
	CONN[NE1][7] = CONN[NE1][4] + nnodesx*nnodesy;
	
	CONN[NE2][1] = NE2;

	CONN[NE2][2] = CONN[NE2 - 2*nelx*nely][5];
	CONN[NE2][3] = CONN[NE2 - 2*nelx*nely][6];
	CONN[NE2][4] = CONN[NE2 - 2*nelx*nely][7];

	CONN[NE2][5] = CONN[NE2][2] + nnodesx*nnodesy;
	CONN[NE2][6] = CONN[NE2][3] + nnodesx*nnodesy;
	CONN[NE2][7] = CONN[NE2][4] + nnodesx*nnodesy;
	
	NE1 = NE2;
      }
    }
  }
	
  if( iblock == 1){
    jstrt = j-1; //Stored from the previous block
  }

  return 0;
}
