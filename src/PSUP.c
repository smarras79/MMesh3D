/******************************************************************
*    *PSUP(int **CONN, int conn_flag, int nnodes, int nelem, int max_nnodes, unsigned int array_numb)
*
*    Points Surrounding Point (PSUP):
*
*    Input:
*    - **CONN is a 2D array containing the mesh. The mesh should have been 
*      already read within the main file that calls PSUP, and should be stored into CONN.
*      Now, CONN can be of two different forms; see below:
*
*    - CONN: array of the connectivity matrix:
*                        WARNING: this function works ONLY for unstructured 2D grids.
*                        It must contain the conn matrix in the format:
*                        conn(1:nelem, 1:elnnode+conn_flag)
*                        where nelem is the mesh number of elements and
*                        elnnode is the max number of nnodes that any element may have.
*                        conn_flag = 1 means that the 1st column of the CONN matrix represents
*                                      the element number (see below);
*                        conn_flag = 0 means that the 1st column of the CONN matrix is already
*                                      a node number (see below);
*                        ex with conn_flag = 1:
*                        1 2 5 3
*                        2 6 3 1
*                        3 3 1 8
*                        4 9 8 5
*                        5 ... ... ...
*                        ...
*                        ex with conn_flag = 0:
*                        2 5 3
*                        6 3 1
*                        3 1 8
*                        9 8 5
*                        ...
*
*
*	- nnodes:     Number of points in the mesh
*	- nelem:      Number of elements in the mesh
*	- max_nnodes: Maximum number of nodes of the elements: ex. if a mesh is made of quadrilaterals,
*                         max_nnodes = 4; if made of triangles, max_nnodes = 3.
*	- array_numb: use 1 of you allocated CONN on a 1-base index (as in Fortran), 
*                  or 0 if you allocated on a 0-index base (standard C)
*    
*    Output (returns):
*    - An array of boundary flags: this array contains the nudes that are on the mountain boundary.
*		WARNING: This array is correct ONLY if the mesh in input is made of quadrilaterals.
*		For this to work on triangles as well, an additional algorithm must be added to this
*		function; the additional algorithm is also found in Lohner [1].
*
*    Revised and fixed: June 20 2009. Working correctly.
*    Simone Marras
*
*    REFERENCES:
*    [1] Lohner,R. (2008), Applied computational fluid dynamics techniques",ed.II, John Wiley and Sons
*
******************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "nrutil.h"

int *PSUP(int **CONN, int conn_flag, int nnodes, int nelem, int max_nnodes, unsigned int array_numb)
{
  
  int bnodes, *BFLAG;
  
  int i,j, ipoin, jpoin, inode, ielem, ipoi1, ipoi2, istor;
  int szesup1;
  int iesup;
  int max_psup;
  
  int *esup1, *esup2, *nesup, *psup1, *psup2, *npsup;
  int *lpoin;
  
  int **CONN_local;
  int TEMP;
  
  //Dynamic allocation
  esup2 = calloc(nnodes+1,sizeof(int*));
  esup1 = calloc(1,sizeof(int*));	// Will be reallocated only later
  nesup = calloc(nnodes+1,sizeof(int*));
  psup2 = calloc(nnodes+1,sizeof(int*));
  lpoin = calloc(nnodes+1,sizeof(int*));
  psup1 = calloc(1,sizeof(int*)); // Will be reallocated only later
  npsup = calloc(nnodes+1,sizeof(int*));
  
  BFLAG = calloc(nnodes+1,sizeof(int *)); // Will be reallocated only later
  //IMPORTANT NOTICE: BFLAG is 0-indexed BUT the valuable indexes for node flagging start with 1!!!
  //					BFLAG[0] has NO MEANING at all.	
  
  if(esup1==NULL || esup2==NULL || nesup==NULL || psup1==NULL || psup2==NULL || npsup==NULL || lpoin==NULL) {
    printf(" # !!! ERROR IN DYNAMIC ALLOCATION\n THE PROGRAM WILL EXIT NOW! ");
    exit(1);
  }
  
  printf("\n###############################\n# Elements surrounding points #\n###############################\n");
  
  /* LOHNER's ALGORITHM BEGINS HERE */
  
  //Element pass 1: count the number of elements connected to each point:
  for(ielem=array_numb; ielem<nelem+array_numb; ielem++){
    for(inode=array_numb; inode<max_nnodes+array_numb; inode++){
      //Update storage counter, storing ahead:
      ipoi1 = CONN[ielem][inode+1] + 1;
      esup2[ipoi1-1+array_numb] = esup2[ipoi1-1+array_numb] + 1;
      //printf("PASS 1 esup2[%d] = %d\n", ipoi1, esup2[ipoi1]); //OK until here
    }
  }
  
  //Storage/reshugffling pass 1:
  for(ipoin=array_numb+1; ipoin<nnodes+1+array_numb; ipoin++){
    //Update storage counter and store:
    esup2[ipoin-1+array_numb] = esup2[ipoin-1+array_numb] + esup2[ipoin-1+array_numb-1]; //Ok
    //printf("esup2[%d] = %d\n", ipoin, esup2[ipoin]); //Ok
  }
  
  szesup1 = esup2[nnodes+array_numb];
  esup1 = realloc(esup1,(szesup1+1)*sizeof(int*));
  
  //Set esup1 to zero:
  for(i=array_numb; i<szesup1+1; i++)
    esup1[i] = 0;
  
  if(esup1 != NULL)
    printf(" #\n # esup1 Now reallocating more memory... \n");
  else
    printf(" !!! NOT ENOUGH MEMORY \n");
  
  //Element pass 2: store the elements in esup1:
  for(ielem=array_numb; ielem<nelem+array_numb; ielem++){
    for(inode=array_numb; inode<max_nnodes+array_numb; inode++){
      //Update storage counter, storing in esup1:
      ipoin = CONN[ielem][inode+1]; //Ok
      //			printf("ipoin = %d\n", ipoin);
      istor = esup2[ipoin-1+array_numb] + 1; //Ok
      //			printf("istor = %d\n", istor);
      
      esup2[ipoin-1+array_numb] = istor+1-array_numb; //OK
      esup1[istor-1+array_numb] = ielem+1-array_numb; //Ok
      
      //			printf("esup2[%d] = %d\n", ipoin-1+array_numb, esup2[ipoin-1+array_numb]); //Ok
      //			printf("esup1[%d] = %d\n", istor-1+array_numb, esup1[istor-1+array_numb]); //Ok
      
      
    }
  }
  
  
  //Storage/reshuffling pass 2 and count the number of elements sorrounding point ipoin:
  for(ipoin = nnodes+array_numb; ipoin>=array_numb+1; ipoin--){
    nesup[ipoin] = esup2[ipoin] - esup2[ipoin-1];
    esup2[ipoin] = esup2[ipoin-1]; //OK
  }
  esup2[array_numb] = 0; //OK
  nesup[array_numb] = esup2[array_numb+1] - esup2[array_numb]; //OK
  
  //for(i=array_numb; i<=nnodes+array_numb; i++){
  //	printf("nesup2[%d] = %d\n", i, nesup[i]);//OK
  //	printf("esup2[%d] = %d\n", i, esup2[i]);//OK
  //}
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Points surrounding points:
  printf("\n ###############################\n # Points surrounding points #\n");
  
  psup2[array_numb] = 0;
  istor = 0;
  
  //!!! IMPORTANT NOTE on the ALGORITHM: Here you will find a set of TWO loops practically identical:
  //the first time that the loop is used, it obtains the number of elements that
  //will fill the array "psup1" that until now was not yet declared and allocated;
  //The second time, it actually fills in such array (this does not happen in the
  //original algorithm of Lohner presented in the book.):
  
  //1) FIRST OF 2 LONG LOOPS to obtain the size of "psup1"
  for(ipoin=array_numb; ipoin<nnodes+array_numb; ipoin++){
    //Loop over the elements sorrounding the point:
    for(iesup=esup2[ipoin]+1; iesup<=esup2[ipoin+1]; iesup++){
      ielem = esup1[iesup-1+array_numb]; //OK
      //printf("iesup = %d  ielem = %d\n", iesup, ielem);//OK
      
      //Loop over the nodes of the element
      for(inode=array_numb; inode<max_nnodes+array_numb; inode++){
	jpoin = CONN[ielem-1+array_numb][inode-1+array_numb+1];
	//printf("jpoin[%d] = %d\n", inode, jpoin); //OK
	
	if(jpoin != ipoin && lpoin[jpoin-1+array_numb] != ipoin){//OK
	  
	  //Update storage counter, storing ahead, and mark lpoin:
	  istor = istor + 1;
	  //printf(" ISTOR = %d\n ", istor);//OK
	  
	  //NOTE: The original loop of the algorithm (corresponding to our 2nd loop below this),
	  //here already fills in psup1! We do it in the next after allocating it.
	  
	  lpoin[jpoin-1+array_numb] = ipoin; //OK
	  //printf(" lpoin[%d] = %d\n ",jpoin-1+array_numb, lpoin[jpoin-1+array_numb]);
	}
      }
    }
  }
  
  
  //Reallocate psup1 with istor elements and BFLAG with n bnodes:
  psup1 = realloc(psup1,(istor+1)*sizeof(int*));		//OK
  lpoin = realloc(lpoin, (nnodes+1)*sizeof(int*));	//OK
  
  if(psup1!=NULL)
    printf(" #\n # psup1 Now reallocating more memory... \n");
  else
    printf("\n !!! NOT ENOUGH MEMORY to allocate psup1\n");
  
  //Null the used arrays before reusing it:
  //OK
  for(i=0; i<nnodes+1; i++){
    lpoin[i] = 0;
    psup2[i] = 0;
    npsup[i] = 0;
  }
  for(i=0; i<istor+1; i++){
    psup1[i] = 0;
  }
  
  istor = 0;
  //2) SECOND OF 2 LONG LOOPS to finally FILL IN "psup1"
  for(ipoin=array_numb; ipoin<nnodes+array_numb; ipoin++){
    //Loop over the elements sorrounding the point:
    for(iesup=esup2[ipoin]+1; iesup<=esup2[ipoin+1]; iesup++){
      ielem = esup1[iesup-1+array_numb]; //OK
      //printf("iesup = %d  ielem = %d\n", iesup, ielem);//OK
      
      //Loop over the nodes of the element
      for(inode=array_numb; inode<max_nnodes+array_numb; inode++){
	jpoin = CONN[ielem-1+array_numb][inode-1+array_numb+1];
	//printf("jpoin[%d] = %d\n", inode, jpoin); //OK
	
	if(jpoin != ipoin && lpoin[jpoin-1+array_numb] != ipoin){//OK
	  
	  //Update storage counter, storing ahead, and mark lpoin:
	  istor = istor + 1;
	  //printf(" ISTOR = %d\n ", istor);//OK
	  psup1[istor-1+array_numb] = jpoin-1+array_numb;//OK
	  //printf(" psup1[%d] = %d\n ",istor-1+array_numb, psup1[istor-1+array_numb]);//OK
	  lpoin[jpoin-1+array_numb] = ipoin-1+array_numb;//OK
	  //printf(" lpoin[%d] = %d\n ",jpoin-1+array_numb, lpoin[jpoin-1+array_numb]);//OK
	  
	}
      }
    }
    
    //Update storage counters:
    psup2[ipoin+1] = istor;//OK
    npsup[ipoin] = psup2[ipoin+1]-psup2[ipoin];//OK
    //printf("psup2[%d] = %d\n", ipoin+1,psup2[ipoin+1]); //OK
    //printf(" # Node %d is surrounded by %d nodes\n", ipoin,npsup[ipoin]); //OK
    
    /*Here the Boundary flag is filled with only the boundary nodes.
      Be very careful because this is correct only for quadrilaterals!!!
      If you want the b.flag to be filled correctly when triangles are used, you
      need to add another part of algorithm following Lohner: chp. 2.2.3, pag.12
      'elements sorrounding elements' 
    */
    
  }//END second loop
  
  //Get the maximum number of points surrounding points:
  // This can be used as the maximum number of nonzero elements
  // appearing among all the rows of the stiffness matrix
  max_psup = 0;
  for(ipoin=0; ipoin<nnodes+1; ipoin++){
    if(npsup[ipoin] > max_psup){
      //printf(" npsup[%d] = %d\n", ipoin, npsup[ipoin]);//OK
      max_psup = npsup[ipoin];//OK
    }
  }
  
  //Here we count the number of nodes that are boundary nodes. This count is
  //necessary to allocate the memory for the Boundary FLAG.
  //IMPORTANT NOTICE: BFLAG is 0-indexed BUT the valuable indexes for node flagging start with 1!!!
  //					BFLAG[0] has NO MEANING at all.
  for(ipoin = 0; ipoin<nnodes+1; ipoin++){
    if( max_nnodes == 3 && npsup[ipoin] < 6)
      //Boundary node for TRI element grids
      {
	BFLAG[ipoin] = 0;
      }
    else if( max_nnodes == 4 &&  npsup[ipoin] < 8)
      //Boundary node for QUAD element grids
      {
	BFLAG[ipoin] = 0;
      }
    else //if internal node
      BFLAG[ipoin] = npsup[ipoin];
    
    //printf(" BFLAG[%d] = %d\n", ipoin, BFLAG[ipoin]);
  }
  
  printf(" #\n # Maximum number of points surrounding points: %d\n###############################\n\n", max_psup);
  
  /*Free memory*/
  
  free(esup2);
  free(esup1);
  free(nesup);
  free(psup2);
  free(lpoin);
  free(psup1);
  free(npsup);
  
  //free(BFLAG);	// NOTICE: Uncomment this if you are not returning BFLAG to the calling function. Otherwise, keep it commented
  // but DO NOT FORGET to free BFLAG in the calling function after use of BFLAG.
  
  return BFLAG;
}
