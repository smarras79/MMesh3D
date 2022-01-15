/*
 * This function builds the connectivity matrix:
 *
 * Simone Marras, Aug 2012
 */
#include <stdlib.h>

#include "MYSTRUCTS.h"
#include "ALMOST_EQUAL.h"
#include "BUILD_CONN.h"
#include "GRID2CONN.h"
#include "GLOBAL_VARS.h"
#include "MEMORY.h"
#include "MYDEFINE.h"
#include "NRUTIL.h"
#include "PRINT.h"
#include "SURFACES.h"

//DEfinitions needed for high order grids
#include "MESH.h"
#include "BUILD_LGL.h"


//BDY CODES for boundary edges:
#define TOP_FLG   111
#define BOTT_FLG  222
#define LEFT_FLG  222
#define RIGHT_FLG 222

#define SOUTH 0
#define NORTH 2
#define EAST  1
#define WEST  3
#define BOTT  4
#define TOP   5

int BUILD_CONN(void)
{
  if( !strncmp(problem[2],"HEXA",4) )
    GRID2CONNhex(0, 0, nnodesx, nnodesy, nnodesz, CONN, ELTYPE, 1);
  else if( !strncmp(problem[2],"WEDGE", 4) )
    GRID2CONNwedge(0, 0, nnodesx, nnodesy, nnodesz, CONN, ELTYPE, 1);
    else{
    printf(" \n");
    printf(" !!! ERROR in input: 'Element_types' can only be HEXA or WEDGE\n" );
    printf(" !!! ERROR in input: HEXA will be set by default\n");
    GRID2CONNhex(0, 0, nnodesx, nnodesy, nnodesz, CONN, ELTYPE, 1);
  }
  
  return 0;
}

/*--------------------------------------------------------------------------
 * Reorder CONN numbering into CGNS:
 *--------------------------------------------------------------------------*/
int CGNS_ORDERING(int **CONN, int nelem)
{
    int CONNcgns[nelem][8];
    
    for (int iel=0; iel<nelem; iel++){	
	CONNcgns[iel][0] = CONN[iel][5];
	CONNcgns[iel][1] = CONN[iel][6];
	CONNcgns[iel][2] = CONN[iel][2];
	CONNcgns[iel][3] = CONN[iel][1];
	CONNcgns[iel][4] = CONN[iel][4];
	CONNcgns[iel][5] = CONN[iel][7];
	CONNcgns[iel][6] = CONN[iel][3];
	CONNcgns[iel][7] = CONN[iel][0];

	//Rewrite:
	CONN[iel][0] = CONNcgns[iel][0];
	CONN[iel][1] = CONNcgns[iel][1];
	CONN[iel][2] = CONNcgns[iel][2];
	CONN[iel][3] = CONNcgns[iel][3];
	CONN[iel][4] = CONNcgns[iel][4];
	CONN[iel][5] = CONNcgns[iel][5];
	CONN[iel][6] = CONNcgns[iel][6];
	CONN[iel][7] = CONNcgns[iel][7];
    }

    return 0;
}

int BUILD_EDGES(int **CONN, int nelem)
{
    /*
     * Get the boundary nodes and the element to which they belong:
     */
    int ibdy_edge   = 0;
    int iface       = 0;
    int iedge_g;
    
    /* Allocate conn_edge_el
     * FACE_LtoG
     * FACE_in_ELEM
     * conn_edge_el
     * conn_face_el
     */
    MEMORY_ALLOCATE(6);

    //VIEW_i2DMAT("CONN", CONN, 0,nelem-1, 0,7);
    
    for (iel = 0; iel<nelem; iel++)
	{
	/*
	 * Edges bottom face:
	 */
	iedge_g = 0;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][0];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][1];
	iedge_g = 1;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][1];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][2];
	iedge_g = 2;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][2];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][3];
	iedge_g = 3;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][3];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][0];
	//Edges top face
	iedge_g = 4;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][4];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][5];
	iedge_g = 5;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][5];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][6];
	iedge_g = 6;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][6];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][7];
	iedge_g = 7;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][7];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][4];
	
	//Vertical edges
	iedge_g = 8;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][0];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][4];
	iedge_g = 9;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][1];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][5];
	iedge_g = 10;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][2];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][6];
	iedge_g = 11;
	conn_edge_el[iel][iedge_g][0] = CONN[iel][3];
	conn_edge_el[iel][iedge_g][1] = CONN[iel][7];
	
	/*
	 * Local faces node connectivity:
	 * i.e. what nodes belong to a given local face in iel:
	 */
	conn_face_el[iel][SOUTH][0] = CONN[iel][0];
	conn_face_el[iel][SOUTH][1] = CONN[iel][1];
	conn_face_el[iel][SOUTH][2] = CONN[iel][5];
	conn_face_el[iel][SOUTH][3] = CONN[iel][4];

	conn_face_el[iel][NORTH][0] = CONN[iel][2];
	conn_face_el[iel][NORTH][1] = CONN[iel][3];
	conn_face_el[iel][NORTH][2] = CONN[iel][7];
	conn_face_el[iel][NORTH][3] = CONN[iel][6];
	
	conn_face_el[iel][EAST][0] = CONN[iel][1];
	conn_face_el[iel][EAST][1] = CONN[iel][2];
	conn_face_el[iel][EAST][2] = CONN[iel][6];
	conn_face_el[iel][EAST][3] = CONN[iel][5];
	
	conn_face_el[iel][WEST][0] = CONN[iel][3];
	conn_face_el[iel][WEST][1] = CONN[iel][0];
	conn_face_el[iel][WEST][2] = CONN[iel][4];
	conn_face_el[iel][WEST][3] = CONN[iel][7];
	
	conn_face_el[iel][BOTT][0] = CONN[iel][3];
	conn_face_el[iel][BOTT][1] = CONN[iel][2];
	conn_face_el[iel][BOTT][2] = CONN[iel][1];
	conn_face_el[iel][BOTT][3] = CONN[iel][0];
	
	conn_face_el[iel][TOP][0] = CONN[iel][4];
	conn_face_el[iel][TOP][1] = CONN[iel][5];
	conn_face_el[iel][TOP][2] = CONN[iel][6];
	conn_face_el[iel][TOP][3] = CONN[iel][7];
	
	/*for (int iface=0; iface<6; iface++)
	    {
		iedge_g = 0;
		conn_face_edge_el[iel][iface][iedge_g][0] = conn_face_el[iel][iface][1];
		conn_face_edge_el[iel][iface][iedge_g][1] = conn_face_el[iel][iface][2];
		conn_face_edge_el[iel][iface][iedge_g][2] = iel;
		
		iedge_g = 1;
		conn_face_edge_el[iel][iface][iedge_g][0] = conn_face_el[iel][iface][2];
		conn_face_edge_el[iel][iface][iedge_g][1] = conn_face_el[iel][iface][3];
		conn_face_edge_el[iel][iface][iedge_g][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge_g][iface][iel][conn_face_edge_el[iel][iface][iedge_g][0]][conn_face_edge_el[iel][iface][iedge_g][1]];
		
		iedge_g = 2;
		conn_face_edge_el[iel][iface][iedge_g][0] = conn_face_el[iel][iface][3];
		conn_face_edge_el[iel][iface][iedge_g][1] = conn_face_el[iel][iface][4];
		conn_face_edge_el[iel][iface][iedge_g][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge_g][iface][iel][conn_face_edge_el[iel][iface][iedge_g][0]][conn_face_edge_el[iel][iface][iedge_g][1]];
		
		iedge_g = 3;
		conn_face_edge_el[iel][iface][iedge_g][0] = conn_face_el[iel][iface][4];
		conn_face_edge_el[iel][iface][iedge_g][1] = conn_face_el[iel][iface][1];
		conn_face_edge_el[iel][iface][iedge_g][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge_g][iface][iel][conn_face_edge_el[iel][iface][iedge_g][0]][conn_face_edge_el[iel][iface][iedge_g][1]];
		}*/
	}

    for(int iel = 0; iel<nelem; iel++) {
	for(int iface = 0; iface<6; iface++) {
	    for(int i = 0; i<2; i++) {
		FACE_in_ELEM[iel][iface][i] = -1;
	    }
	}
    }
    
    /*
     * sort by 3rd dimension
     */
    int ***conn_face_el_sort;
    int auxi[5];
    int npoints_in_iface = 4;
    conn_face_el_sort = i3tensor(0, nelem-1, 0, 5, 0, 3);
    for (int iel=0; iel<nelem; iel++) {
	for (int j=0; j<6; j++) {
	    for (int k=0; k<NELFACES; k++) {
		int kk = k + 1;
		auxi[kk] = conn_face_el[iel][j][k];
	    }	    
	    isort(npoints_in_iface, auxi);
	    for (int k=0; k<4; k++) {
		int kk = k + 1;
		conn_face_el_sort[iel][j][k] = auxi[kk];
	    }
	}
    }
    //VIEW_i3DMAT("CONN_FACE_EL_SORT", conn_face_el_sort, 0, nelem-1, 0, 5, 0, 3);
    //printf(" END SORT\n");
    //END sorting
    
    iface = 0;
    for(int iel=0; iel<nelem; iel++) {
	for(int jel=0; jel<nelem; jel++) {

	    if(     conn_face_el_sort[iel][BOTT][0] == conn_face_el_sort[jel][TOP][0] && \
		    conn_face_el_sort[iel][BOTT][1] == conn_face_el_sort[jel][TOP][1] && \
		    conn_face_el_sort[iel][BOTT][2] == conn_face_el_sort[jel][TOP][2] && \
		    conn_face_el_sort[iel][BOTT][3] == conn_face_el_sort[jel][TOP][3])
		{	
		    FACE_in_ELEM[iel][BOTT][0] = iel;
		    FACE_in_ELEM[iel][BOTT][1] = jel;
				
		    FACE_in_ELEM[jel][TOP][0] = jel;
		    FACE_in_ELEM[jel][TOP][1] = iel;

		    //printf(" SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n", BOTT+1,iel+1, TOP+1, jel+1);
						
		    iface = iface + 1;
			
		} else if (conn_face_el_sort[iel][EAST][0] == conn_face_el_sort[jel][WEST][0] && \
			   conn_face_el_sort[iel][EAST][1] == conn_face_el_sort[jel][WEST][1] && \
			   conn_face_el_sort[iel][EAST][2] == conn_face_el_sort[jel][WEST][2] && \
			   conn_face_el_sort[iel][EAST][3] == conn_face_el_sort[jel][WEST][3])
		{
			
		    FACE_in_ELEM[iel][EAST][0] = iel;
		    FACE_in_ELEM[iel][EAST][1] = jel;
			
		    FACE_in_ELEM[jel][WEST][0] = jel;
		    FACE_in_ELEM[jel][WEST][1] = iel;
			
		    //printf(" SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n", EAST+1, iel+1, WEST+1, jel+1);
				
		    iface = iface + 1;
			
		} else if  (conn_face_el_sort[iel][SOUTH][0] == conn_face_el_sort[jel][NORTH][0] && \
			    conn_face_el_sort[iel][SOUTH][1] == conn_face_el_sort[jel][NORTH][1] && \
			    conn_face_el_sort[iel][SOUTH][2] == conn_face_el_sort[jel][NORTH][2] && \
			    conn_face_el_sort[iel][SOUTH][3] == conn_face_el_sort[jel][NORTH][3])
		{
			
		    FACE_in_ELEM[iel][SOUTH][0] = iel;
		    FACE_in_ELEM[iel][SOUTH][1] = jel;
			
		    FACE_in_ELEM[jel][NORTH][0] = jel;
		    FACE_in_ELEM[jel][NORTH][1] = iel;
			
		    //printf(" SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n", SOUTH+1, iel+1, NORTH+1, jel+1);

		    iface = iface + 1;
		}
	}
    }
    int face_direction[6] = {SOUTH, EAST, NORTH, WEST, BOTT, TOP};
    int nint_faces        = iface;
    nfaces                = nint_faces + nbdy_faces;
    
    /*----------------------------------------------------------------------------------------------------
     * Store connectivity of internal faces: 
     * start from the boundary faces and append the internal faces afterwards:
     *----------------------------------------------------------------------------------------------------*/
    //1) Store bdy faces in CONN_FACE(1:BDY_FACES, 1:4)
    //int CONN_FACE[nfaces][4];
    int nface_all = nelem*6;
    int CONN_FACE_all[nface_all][4];

    MEMORY_ALLOCATE(7);
    
    //printf(" CONN_BDY_FACES\n");
    for (int i=0; i<nbdy_faces; i++) {
	for (int j=0; j<4; j++) {
	    CONN_FACE[i][j] = CONN_BDY_FACES[i][j];
	    //printf("%d ", CONN_FACE[i][j]);
	}
	//printf("\n");
    }
    
    /*
     * 2) Append the internal faces to CONN_FACE(NCONN_BDY_FACES+1:NFACES, 1:4)
     *
     *
     * Loop 1 to populate internal edges to CONN_EDGE:
     *
     *
     * Detect ALL faces, including the shared ones:
     */
    int iface_all = 0;
    for (int iel=0; iel<nelem; iel++) {
	for (int ifac=0; ifac<6; ifac++) {
	    CONN_FACE_all[iface_all][0] = conn_face_el[iel][ifac][0];
	    CONN_FACE_all[iface_all][1] = conn_face_el[iel][ifac][1];
	    CONN_FACE_all[iface_all][2] = conn_face_el[iel][ifac][2];
	    CONN_FACE_all[iface_all][3] = conn_face_el[iel][ifac][3];
	    
	    iface_all = iface_all + 1;
	}
    }
    
    int CONN_FACE_tmp[iface_all][4];
    for (int i=0; i<iface_all; i++) {
	for (int j=0; j<4; j++) {
	    int jj = j + 1;
	    auxi[jj] = CONN_FACE_all[i][j];
	}
	isort(4, auxi);
	for (int j=0; j<4; j++) {
	    int jj = j + 1;
	    CONN_FACE_tmp[i][j] = auxi[jj];
	}
    }//OK
    
    /*for (int i=0; i<iface_all; i++) {
	printf(" %d %d %d %d %d\n", i, CONN_FACE_tmp[i][0], CONN_FACE_tmp[i][1], CONN_FACE_tmp[i][2],  CONN_FACE_tmp[i][3]);
	}*/
    //int nface_all = iface_all;  /* NOTICE: nface_all are ALL of the faces,
    //				     * including the repeated faces. 
    //				     * nfaces are only the faces counted once.*/
    int REPEATED_auxi[nfaces][2];
    int REPEATED_index[nface_all+1];
    int FACE_MULTIPLICITY_auxi[nface_all+1];
    for (int i=0; i<nfaces; i++) {
	REPEATED_auxi[i][0] = -1;
	REPEATED_auxi[i][1] = -1;
    }    
    for (int i=0; i<nface_all; i++) {
	REPEATED_index[i] = -1;
	FACE_MULTIPLICITY_auxi[i] = 0;
    }
    int iface_repeated = 0;
    int krepeated      = 0;

    printf(" # Number of all faces\t\t%d\n # number of unique faces\t%d\n # Number of interlan faces\t%d\n # Number of boundary faces\t%d\n", nface_all, nfaces, nint_faces, nbdy_faces);
    
    iface = 0;
    for (int i=1; i<=nface_all; i++) {
	int multiplicity = 0;
	//	printf("Face %d -> [%d %d %d %d] \n", i, CONN_FACE_tmp[i-1][0], CONN_FACE_tmp[i-1][1], CONN_FACE_tmp[i-1][2], CONN_FACE_tmp[i-1][3]);	
	for (int j=i; j<=nface_all; j++) {
	    if (j != i) {
		if ( iAlmostEqual(CONN_FACE_tmp[i-1][0], CONN_FACE_tmp[j-1][0]) && \
		     iAlmostEqual(CONN_FACE_tmp[i-1][1], CONN_FACE_tmp[j-1][1]) && \
		     iAlmostEqual(CONN_FACE_tmp[i-1][2], CONN_FACE_tmp[j-1][2]) && \
		     iAlmostEqual(CONN_FACE_tmp[i-1][3], CONN_FACE_tmp[j-1][3]) ) {
		    
		    multiplicity = multiplicity  + 1;
		    krepeated = krepeated + 1;
		    REPEATED_index[krepeated-1] = j;
		    
		    FACE_MULTIPLICITY_auxi[nbdy_faces + iface] = multiplicity + 1;

		    CONN_FACE[nbdy_faces + iface][0] =  CONN_FACE_tmp[j-1][0];
		    CONN_FACE[nbdy_faces + iface][1] =  CONN_FACE_tmp[j-1][1];
		    CONN_FACE[nbdy_faces + iface][2] =  CONN_FACE_tmp[j-1][3];
		    CONN_FACE[nbdy_faces + iface][3] =  CONN_FACE_tmp[j-1][2];
		    
		    /*printf("k=%d, REPEATED_index %d -> [%d %d %d %d] repeated %d times (%d)\n", krepeated-1, REPEATED_index[krepeated-1], \
			   CONN_FACE_tmp[j-1][0], CONN_FACE_tmp[j-1][1], CONN_FACE_tmp[j-1][2], CONN_FACE_tmp[j-1][3], \
			   FACE_MULTIPLICITY_auxi[nbdy_faces + iface], multiplicity + 1);
		    */
		    iface = iface + 1;
		}
	    }
	}
    } //OK (DO NOT CALL "CGNS_ORDERING(CONN, nelem)"!!!!!)
    //printf(" NBDY FACEs + IFACE = %d\n", nbdy_faces + iface);
    //VIEW_i2DMAT("CONN_FACE", CONN_FACE, 0, nfaces-1, 0, 3);
    
    /*
     * Set repeated entries of CONN_FACE_all to -1:
     */
    int nrepeated = krepeated;
    for (int i=0; i<nrepeated; i++) {
	int irepeated_index = REPEATED_index[i];
	CONN_FACE_all[irepeated_index][0] = -1;
	CONN_FACE_all[irepeated_index][1] = -1;
	CONN_FACE_all[irepeated_index][2] = -1;
	CONN_FACE_all[irepeated_index][3] = -1;
    }
    
    /* /\* THIS WAS THE OLD WAY AS IN THE MATLAB CODE BUT IT'S NOT NECESSARY.*\/ */
    /* /\*  * STORE EACH FACE into CONN_FACE(1:`, 1:4): *\/ */
    /* /\*  * *\/ */
    /* /\* int FACE_MULTIPLICITY[nface_all]; //overallocated *\/ */
    /* /\* int ifac = 0; *\/ */
    /* /\* for (int iface_all=0; iface_all<nface_all; iface_all++) { *\/ */
    /* /\* 	//if (CONN_FACE_all[iface_all][0] > 0) { *\/ */
    /* /\* 	if (FACE_MULTIPLICITY_auxi[nbdy_faces+iface_all] > 0) { *\/ */
		    
    /* /\* 	    //CONN_FACE[ifac][0] = CONN_FACE_all[iface_all][0]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][1] = CONN_FACE_all[iface_all][1]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][2] = CONN_FACE_all[iface_all][2]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][3] = CONN_FACE_all[iface_all][3]; *\/ */

    /* /\* 	    //CONN_FACE[ifac][0] = CONN_FACE_tmp[iface_all][0]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][1] = CONN_FACE_tmp[iface_all][1]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][2] = CONN_FACE_tmp[iface_all][2]; *\/ */
    /* /\* 	    //CONN_FACE[ifac][3] = CONN_FACE_tmp[iface_all][3]; *\/ */
	    
    /* /\* 	    FACE_MULTIPLICITY[nbdy_faces + ifac] = FACE_MULTIPLICITY_auxi[iface_all]; *\/ */
	
    /* /\* 	    // FACE_in_ELEM[ifac][0] = FACE_in_ELEM[iface_all][0]; *\/ */
    /* /\* 	    // FACE_in_ELEM[ifac][1] = FACE_in_ELEM[iface_all][1]; *\/ */
	    
    /* /\* 	    ifac = ifac + 1; *\/ */
    /* /\* 	} *\/ */
    /* /\* 	}*\\/ *\/ */


    //OK WORKING: CONN_FACE is now correctly populated and ordered.
    //Uncomment the following for loop to print it to screen.
    /*for (int iface=0; iface<nfaces; iface++) {
	if (iface < nbdy_faces) {
	    printf(" BDY: CONN_FACE(%d,1:4) = %d %d %d %d - repeated %d times\n", iface, CONN_FACE[iface][0], CONN_FACE[iface][1], CONN_FACE[iface][2], CONN_FACE[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	} else {
	    printf(" INT: CONN_FACE(%d,1:4) = %d %d %d %d - repeated %d times\n", iface, CONN_FACE[iface][0], CONN_FACE[iface][1], CONN_FACE[iface][2], CONN_FACE[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	}
	} //OK CONN_FACE and FACE_MULTIPLICITY;*/
    
    
    /*--------------------------------------------------------------------------
     * Populate FACE_LtoG(1:NEL,1:6) OK
     *--------------------------------------------------------------------------*/
    int CONN_FACE_sort[nfaces][4];
    for (int i=0; i<nfaces; i++) {
	for (int j=0; j<4; j++) {
	    int jj = j + 1;
	    auxi[jj] = CONN_FACE[i][j];  //THIS IS PROBABLY INCORRECTLY DONE! RECHECK IT
	}
	isort(npoints_in_iface, auxi);
	for (int j=0; j<4; j++) {
	    int jj = j + 1;
	    CONN_FACE_sort[i][j] = auxi[jj];
	}
    }//OK

    /*for (int iface=0; iface<nfaces; iface++) {
	if (iface < nbdy_faces) {
	    printf(" face %d: s BDY: CONN_FACE_sort(%d,1:4) = %d %d %d %d - repeated %d times\n", iface, iface, CONN_FACE_sort[iface][0], CONN_FACE_sort[iface][1], CONN_FACE_sort[iface][2], CONN_FACE_sort[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	    printf(" face %d: u BDY: CONN_FACE     (%d,1:4) = %d %d %d %d - repeated %d times\n", iface, iface, CONN_FACE[iface][0], CONN_FACE[iface][1], CONN_FACE[iface][2], CONN_FACE[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	} else {
	    printf(" face %d: s INT: CONN_FACE_sort(%d,1:4) = %d %d %d %d - repeated %d times\n", iface, iface, CONN_FACE_sort[iface][0], CONN_FACE_sort[iface][1], CONN_FACE_sort[iface][2], CONN_FACE_sort[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	    printf(" face %d: u INT: CONN_FACE     (%d,1:4) = %d %d %d %d - repeated %d times\n", iface, iface, CONN_FACE[iface][0], CONN_FACE[iface][1], CONN_FACE[iface][2], CONN_FACE[iface][3], FACE_MULTIPLICITY_auxi[iface]);
	}
	}*/ //OK
    
    int IBDY_FACE = 0;
    for (int IFACE=0; IFACE<nfaces; IFACE++) {
	for (int iel=0; iel<nelem; iel++) {
	    for (int iface=0; iface<6; iface++) {
		
		if ( (  CONN_FACE_sort[IFACE][0] == conn_face_el_sort[iel][iface][0] && \
			CONN_FACE_sort[IFACE][1] == conn_face_el_sort[iel][iface][1] && \
			CONN_FACE_sort[IFACE][2] == conn_face_el_sort[iel][iface][3] && \
			CONN_FACE_sort[IFACE][3] == conn_face_el_sort[iel][iface][2]) ||
		      (  CONN_FACE_sort[IFACE][0] == conn_face_el_sort[iel][iface][0] && \
			 CONN_FACE_sort[IFACE][1] == conn_face_el_sort[iel][iface][1] && \
			 CONN_FACE_sort[IFACE][2] == conn_face_el_sort[iel][iface][2] && \
			 CONN_FACE_sort[IFACE][3] == conn_face_el_sort[iel][iface][3])
		     ) {
			    
		    FACE_LtoG[iel][iface] = IFACE;
		    /*printf("  --- FACE_LtoG[%d,%d] = %d -> [%d %d %d %d] \n", iel+1, iface+1, \
			   FACE_LtoG[iel][iface], conn_face_el_sort[iel][iface][0], \
			   conn_face_el_sort[iel][iface][1],		\
			   conn_face_el_sort[iel][iface][2],		\
			   conn_face_el_sort[iel][iface][3]);*/
		}
	    }
	}

	/* UNCOMMENT this part to populate CONN_BDY_FACE[0:NBDY_FACES][0:3]
	   if it is not coming from the GMSH file 
	   
	   if (FACE_MULTIPLICITY_auxi[IFACE] == 0) { //CHECK THE SIZE OF THIS FIRST
	      CONN_BDY_FACE[IBDY_FACE][0] = CONN_FACE[IFACE][0];
	      CONN_BDY_FACE[IBDY_FACE][1] = CONN_FACE[IFACE][1];
	      CONN_BDY_FACE[IBDY_FACE][2] = CONN_FACE[IFACE][2];
	      CONN_BDY_FACE[IBDY_FACE][3] = CONN_FACE[IFACE][3];
	      IBDY_FACE = IBDY_FACE + 1;
           }
	*/
    }
    free_i3tensor(conn_face_el_sort, 0, nelem-1, 0, 5, 0, 3);

    /*----------------------------------------------------------------------------------------------------------------
     * Global edge connectivity:
     *----------------------------------------------------------------------------------------------------------------
     *
     * Loop 1 to populate internal edges to CONN_EDGE:
     *

     *EDGE_FLG(1:NEDGES)          = -1;
     *EDGE_LtoG(1:NELEM, 1:12)    = -1; %map from local edge number (1:4) to global edge number
     *EDGE_in_ELEM(1:NEDGES, 1:2) = -1;
     *CONN_EDGE(1:NEDGES, 1:3)    = -1; %third column is to color the edge if repeated.
     *EDGE_MULTIPLICITY(1:NEDGES) = -1;
     *
     * Detect ALL edges, including the shared ones:
     */    
    int nedge_all = nelem*12;
    int CONN_EDGE_tmp[nedge_all][2];
    int CONN_EDGE_all[nedge_all][2];
    int EDGE_in_ELEM[ nedge_all];
    int tmp1;
    int iedge_all = 0;
    for (int iel=0; iel<nelem; iel++) {
	for (int iedg=0; iedg<12; iedg++) {

	     /*
	     * Order first and second point in ascending order [for comparison
	     * purposes only)
	     */
	    if (conn_edge_el[iel][iedg][0] > conn_edge_el[iel][iedg][1]) {
		tmp1 = conn_edge_el[iel][iedg][1];
		conn_edge_el[iel][iedg][1] = conn_edge_el[iel][iedg][0];
		conn_edge_el[iel][iedg][0] = tmp1;
	    }
	    
	    CONN_EDGE_all[iedge_all][0] = conn_edge_el[iel][iedg][0];
	    CONN_EDGE_all[iedge_all][1] = conn_edge_el[iel][iedg][1];
		
	    //Used later to detect repeated rows only:
	    CONN_EDGE_tmp[iedge_all][0] = conn_edge_el[iel][iedg][0];
	    CONN_EDGE_tmp[iedge_all][1] = conn_edge_el[iel][iedg][1];
	    
	    EDGE_in_ELEM[iedge_all] = iel;
	    
	    //printf("iedge = %d (%d, %d) is in elements %d\n", iedge_all+1,  CONN_EDGE_tmp[iedge_all][0], CONN_EDGE_tmp[iedge_all][1], EDGE_in_ELEM[iedge_all]);
	    
	    iedge_all = iedge_all + 1;
	}
    }
    
    krepeated                        =  0;
    int IEDGE_repeated, multiplicity =  0;
    int EDGE_MULTIPLICITY_auxi[nedge_all];
    int EDGE_MULTIPLICITY[nedge_all];      //This is over-allocated but it will do.
    int EDGE_REPEATED_index[nedge_all];
    
    for (int i=1; i<=nedge_all; i++) {
	EDGE_MULTIPLICITY_auxi[i] = 0;
	EDGE_MULTIPLICITY[i] = 0;
    }
    
    for (int i=1; i<=nedge_all; i++) {
	multiplicity = 0;
	for (int j=i; j<=nedge_all; j++) {
	    
	    if (j != i) {
		if ( iAlmostEqual(CONN_EDGE_tmp[i-1][0], CONN_EDGE_tmp[j-1][0]) && \
		     iAlmostEqual(CONN_EDGE_tmp[i-1][1], CONN_EDGE_tmp[j-1][1]) ) {
		    
		    multiplicity                       = multiplicity + 1;
		    krepeated                          = krepeated + 1;
		    EDGE_REPEATED_index[krepeated - 1] = j;
		    EDGE_MULTIPLICITY_auxi[i - 1]      = multiplicity + 1;

		    //printf("i = %d (%d, %d) - repeated %d times\n", i, CONN_EDGE_tmp[i-1][0] , CONN_EDGE_tmp[i-1][1], EDGE_MULTIPLICITY_auxi[i-1]);
		}	
	    }
	}
    } //OK
    
    /*
     * Set repeated entries of CONN_EDGE_tmp to -1:
     */
    int irepeated_index;
    nrepeated = krepeated;
    for (int i=0; i<nrepeated; i++) {
	
	irepeated_index = EDGE_REPEATED_index[i];
	//printf("%d irepeated_index = %d\n", i+1, irepeated_index);
	CONN_EDGE_all[irepeated_index-1][0] = -1;
	CONN_EDGE_all[irepeated_index-1][1] = -1;
    }
    
    /*--------------------------------------------------------------------------
     * STORE EACH EDGE into CONN_EDGE(1:NEDGES, 1:2): OK
     *--------------------------------------------------------------------------*/
    //Count unique edges
    nedges = 0;
    for (int iedge_all=0; iedge_all<nedge_all; iedge_all++){
	if (CONN_EDGE_all[iedge_all][0] > 0) {
	    nedges++;
	}
    }    
    
    //if (irank == ) printf(" N. unique edges (nedges) = %d\n", nedges);
    printf(" # Number of unique edges\t%d\n", nedges);
    
    MEMORY_ALLOCATE(8); //allocate CONN_EDGES[nedges][nop+1];
        
    iedge_g = 0;
    iedge_all = 0;
    for (int iel=0; iel<nelem; iel++) {
	for (int iedg=0; iedg<12; iedg++) {
	    
	    if (CONN_EDGE_all[iedge_all][0] > 0) {		
		CONN_EDGE[iedge_g][0]      = CONN_EDGE_all[iedge_all][0];
		CONN_EDGE[iedge_g][1]      = CONN_EDGE_all[iedge_all][1];
		EDGE_MULTIPLICITY[iedge_g] = EDGE_MULTIPLICITY_auxi[iedge_all];
		//printf("  IEDGE %d has multiplicity %d = (%d, %d)\n", iedge_g+1, EDGE_MULTIPLICITY[iedge_g],CONN_EDGE[iedge_g][0],CONN_EDGE[iedge_g][1]);
		iedge_g++;
	    }
	    iedge_all++;
	}
    }//OK EDGE_MULTIPLICITY[iedge_g], CONN_EDGE[iedge_g][0,1] are correctly populated
   
    /*--------------------------------------------------------------------------
     * Populate EDGE_LtoG(1:nelem,1:12) OK
     *--------------------------------------------------------------------------*/
    for (int iedge_g=0; iedge_g<nedges; iedge_g++){
	for (int iel=0; iel<nelem; iel++) {
	    for (int iedg_el=0; iedg_el<12; iedg_el++) {
		
		if ( (CONN_EDGE[iedge_g][0] == conn_edge_el[iel][iedg_el][0] && CONN_EDGE[iedge_g][1] == conn_edge_el[iel][iedg_el][1]) || \
		     (CONN_EDGE[iedge_g][1] == conn_edge_el[iel][iedg_el][0] && CONN_EDGE[iedge_g][0] == conn_edge_el[iel][iedg_el][1]) ) {
		    
		    EDGE_LtoG[iel][iedg_el] = iedge_g;
		    //printf("  --- EDGE_LtoG(iel=%d,iedg=%d) = %d with nodes --> (%d, %d) with multiplicity %d\n", iel+1, iedg_el+1, EDGE_LtoG[iel][iedg_el],  CONN_EDGE[iedge_g][0], CONN_EDGE[iedge_g][1], EDGE_MULTIPLICITY[iedge_g]);
		}
	    }
	}
    }//OK

    
    /*-------------------------------------------------------------------------
     * BUILD BDY_EDGE(1:NBDY_EDGES, 1:5): 
     * NOTICE: This way of identifying boundary edges only works for structured
     * grids with hexas.
     *-------------------------------------------------------------------------*
    int IBDY_EDGE = 0;
    for (int iedge_g=0; iedge_g<nedges; iedge_g++){
	if (EDGE_MULTIPLICITY[iedge_g] <= 3) {
	    
		 IBDY_EDGE = IBDY_EDGE + 1;
		 
		 BDY_EDGE[IBDY_EDGE][0] = iedge_g;
		 BDY_EDGE[IBDY_EDGE][1] = CONN_EDGE[iedge_g][0];
		 BDY_EDGE[IBDY_EDGE][2] = CONN_EDGE[iedge_g][1];
		 BDY_EDGE[IBDY_EDGE][3] = -1; //IEL;
		 BDY_EDGE[IBDY_EDGE][4] = 1;
		 
		 printf("  edge %d is a BDY EDGE with nodes (%d %d)\n", IEDGE, BDY_EDGE[IBDY_EDGE][1], BDY_EDGE[IBDY_EDGE][2]);
	}
    }
    /*-------------------------------------------------------------------------
     * END BUILD BDY_EDGE(1:NBDY_EDGES, 1:5)
     *-------------------------------------------------------------------------*/
        
    return 0;
}


int ADD_HIGH_ORDER_NODES(void)
{
    if (nop < 2) return 0;
    
    printf(" #------------------------------------------------------------------#\n");
    printf(" # POPULATE GRID with SPECTRAL NODES............................\n");
    
    /*--------------------------------------------------------------------------
     * Auxiliary elemental parameters/counters:
     *--------------------------------------------------------------------------*/
    st_lgl lgl;
    
    int nnodes_linear            = nnodes;
    int ncorner_nodes            = 8;
    int nvolume_internal_nodes   = nelem*(ngl-2)*(ngl-2)*(ngl-2);
    int tot_edges_internal_nodes = nedges*(ngl-2);
    int tot_faces_internal_nodes = nfaces*(ngl-2)*(ngl-2);
    int iedge_g = 1;

    //Update nnodes to total nnodes for high order grids
    nnodes = nnodes_linear + tot_edges_internal_nodes + tot_faces_internal_nodes + nvolume_internal_nodes;

    printf(" # Total number of unique linear nodes    %d\n", nnodes_linear);
    printf(" # Total number of unique nodes on edges  %d\n", tot_edges_internal_nodes);
    printf(" # Total number of unique nodes on faces  %d\n", tot_faces_internal_nodes);
    printf(" # Total number of unique nodes in volume %d\n", nvolume_internal_nodes);
    printf(" # Total number of unique high-order nodes %d\n", nnodes);
    
    //Create coordinates and weights of LGL points along 1D reference element of order nop. E.g. o-x--x-o for nop=3
    lgl = BUILD_LGL(nop);
    
    FILE *fileid, *fileidHO_edges,  *fileidHO_faces, *fileidHO_vol;
    
    fileid         = fopen("COORDS_LO.dat", "w");
    fileidHO_edges = fopen("COORDS_HO_edges.dat", "w");
    fileidHO_faces = fopen("COORDS_HO_faces.dat", "w");
    fileidHO_vol   = fopen("COORDS_HO_vol.dat", "w");

    /*
     * 1. Allocate COORDS_HO
     * 2. store  COORDS_HO[1:nnodes_linear] <- COORDS[1:nnodes_linear]
     * 3. Reeallocate COORDS[1:nnodes] where nnodes is the new nnodes for the HO grid
     */
    MEMORY_ALLOCATE(11);
    for(int ip = 0; ip<nnodes_linear; ip++) {
	COORDS_HO[ip][0] = COORDS[ip][0];
	COORDS_HO[ip][1] = COORDS[ip][1];
	COORDS_HO[ip][2] = COORDS[ip][2];
	fprintf(fileid, " %f %f %f\n", COORDS_HO[ip][0], COORDS_HO[ip][1], COORDS_HO[ip][2]);
    }
    
    /*--------------------------------------------------------------------------
     * Build high order grid points on every edges
     --------------------------------------------------------------------------*/
    double x1, y1, z1;
    double x2, y2, z2;
    double xi;
    int ip1, ip2;
    
    int ip = nnodes_linear + 1; //we start populating from the low order numbering
    for(int iedge_g = 0; iedge_g<nedges; iedge_g++) {
	
	ip1 = min(CONN_EDGE[iedge_g][0], CONN_EDGE[iedge_g][1]) - 1;
	ip2 = max(CONN_EDGE[iedge_g][0], CONN_EDGE[iedge_g][1]) - 1;
	
	x1  = COORDS[ip1][0], y1  = COORDS[ip1][1], z1  = COORDS[ip1][2];
	x2  = COORDS[ip2][0], y2  = COORDS[ip2][1], z2  = COORDS[ip2][2];
	//printf(" IP1 %d = (%f %f %f) --- IP2 %d = (%f %f %f)\n", ip1+1,  x1, y1, z1, ip2+1, x2, y2, z2);
	//OK
	
	/*
	 * Plot LGL points on edges:
	 */
	for(int l=1; l<ngl-1; l++) {
	    
	    //printf(" #\t LGL nodes: X_LG[%d] = %.16f; w_LG[%d] = %.16f\n", l, lgl.ksi[l], l, lgl.weights[l]);
	    
	    xi = lgl.ksi[l];
	    
	    COORDS_HO[ip][0] = x1*(1.0 - xi)*0.5 + (1.0 + xi)*x2*0.5;
	    COORDS_HO[ip][1] = y1*(1.0 - xi)*0.5 + (1.0 + xi)*y2*0.5;
	    COORDS_HO[ip][2] = z1*(1.0 - xi)*0.5 + (1.0 + xi)*z2*0.5;
	    
	    fprintf(fileidHO_edges, " %d %f %f %f\n", ip, COORDS_HO[ip][0], COORDS_HO[ip][1], COORDS_HO[ip][2]);
	    //printf(" ngl %d, iedge %d, ilgl %d, %d %f %f %f\n", ngl, iedge_g, l, ip, COORDS_HO[ip][0], COORDS_HO[ip][1], COORDS_HO[ip][2]);

	    //iconn = iconn + 1;	    
	    ip = ip + 1; //Initialized to highest low order value of npoin.
	}
    }
    
    /*--------------------------------------------------------------------------
     * Populate Faces with high-order points:
     *--------------------------------------------------------------------------*/
    double xa, ya, za;
    double xb, yb, zb;
    double xc, yc, zc;
    double xd, yd, zd;
    double zeta;
        
    int iconn, iconn_face_internal;
    int ip3, ip4;

    ip = nnodes_linear + tot_edges_internal_nodes + 1;
    for(int iface=0; iface<nfaces; iface++) {
	    	
	/*--------------------------------------------------------------------------
	 * Corners ordered counter-clockwise
	 *--------------------------------------------------------------------------*/
	ip1 = CONN_FACE[iface][0]-1;
	ip2 = CONN_FACE[iface][1]-1;
	ip3 = CONN_FACE[iface][2]-1;
	ip4 = CONN_FACE[iface][3]-1;
	//printf(" ip1, ip2, ip3, ip4 = %d %d %d %d\n", ip1+1, ip2+1, ip3+1, ip4+1);
	
	xa = COORDS[ip1][0]; ya = COORDS[ip1][1]; za = COORDS[ip1][2];
	xb = COORDS[ip2][0]; yb = COORDS[ip2][1]; zb = COORDS[ip2][2];
	xc = COORDS[ip3][0]; yc = COORDS[ip3][1]; zc = COORDS[ip3][2];
	xd = COORDS[ip4][0]; yd = COORDS[ip4][1]; zd = COORDS[ip4][2];
	
	//iconn = iconnCurrent + 1;
	//iconn_face_internal = 0;
	for(int k=1; k<ngl-1; k++) {
	    zeta = lgl.ksi[k];

	    for(int i=1; i<ngl-1; i++) {
		xi = lgl.ksi[i];
		
		COORDS_HO[ip][0] = xa*(1 - xi)*(1 - zeta)*0.25 +	\
		    xb*(1 + xi)*(1 - zeta)*0.25 +			\
		    xc*(1 + xi)*(1 + zeta)*0.25 +			\
		    xd*(1 - xi)*(1 + zeta)*0.25;
		
		
		COORDS_HO[ip][1] = ya*(1 - xi)*(1 - zeta)*0.25 +	\
		    yb*(1 + xi)*(1 - zeta)*0.25 +			\
		    yc*(1 + xi)*(1 + zeta)*0.25 +			\
		    yd*(1 - xi)*(1 + zeta)*0.25;
		
		
		COORDS_HO[ip][2] = za*(1 - xi)*(1 - zeta)*0.25 +	\
		    zb*(1 + xi)*(1 - zeta)*0.25 +			\
		    zc*(1 + xi)*(1 + zeta)*0.25 +			\
		    zd*(1 - xi)*(1 + zeta)*0.25;
		
		fprintf(fileidHO_faces, " %d %f %f %f\n", ip, COORDS_HO[ip][0], COORDS_HO[ip][1], COORDS_HO[ip][2]);
			    
		/*
		 * Add internal LGL points to CONN:
		 */
		//FACE_POINT_CONN[iface][i][k] = IP;
				    
		ip = ip + 1;
		//iconn = iconn + 1;
		//iconn_face_internal = iconn_face_internal + 1;
	    }
	}
    }

    /*--------------------------------------------------------------------------
     * Populate internal/Volume high-order points:
     *--------------------------------------------------------------------------*/
    int ip5, ip6, ip7, ip8;
    
    double xe, ye, ze;
    double xf, yf, zf;
    double xg, yg, zg;
    double xh, yh, zh;
    double eta;
    
    //int ivol_start         = nnodes_linear + nedge_internal_nodes*12 + nface_internal_nodes*6 + 1;
    //int iconn_volume_index = nnodes_linear + nedge_internal_nodes*12 + nface_internal_nodes*6 + 1;
    
    ip = nnodes_linear + tot_edges_internal_nodes + tot_faces_internal_nodes + 1;
    ip = 0;
    for (int iel=0; iel<nelem; iel++) {
	
	/*printf(" CONN[%d] = (%d %d %d %d %d %d %d %d)\n", iel, CONN[iel][0],CONN[iel][1],CONN[iel][2],CONN[iel][3],CONN[iel][4], \
	  CONN[iel][5],CONN[iel][6],CONN[iel][7]);*/
	
	/*--------------------------------------------------------------------------
	 * Volume:
	 *--------------------------------------------------------------------------*/
	//iconn  = iconn_volume_index + 1;

	ip1 = CONN[iel][0]-1; ip2 = CONN[iel][1]-1; ip3 = CONN[iel][2]-1; ip4 = CONN[iel][3]-1;
	ip5 = CONN[iel][4]-1; ip6 = CONN[iel][5]-1; ip7 = CONN[iel][6]-1; ip8 = CONN[iel][7]-1;
	
	xa = COORDS[ip1][0]; ya = COORDS[ip1][1]; za = COORDS[ip1][2];
	xb = COORDS[ip2][0]; yb = COORDS[ip2][1]; zb = COORDS[ip2][2];
	xc = COORDS[ip3][0]; yc = COORDS[ip3][1]; zc = COORDS[ip3][2];
	xd = COORDS[ip4][0]; yd = COORDS[ip4][1]; zd = COORDS[ip4][2];
	xe = COORDS[ip5][0]; ye = COORDS[ip5][1]; ze = COORDS[ip5][2];
	xf = COORDS[ip6][0]; yf = COORDS[ip6][1]; zf = COORDS[ip6][2];
	xg = COORDS[ip7][0]; yg = COORDS[ip7][1]; zg = COORDS[ip7][2];
	xh = COORDS[ip8][0]; yh = COORDS[ip8][1]; zh = COORDS[ip8][2];

	iconn = ncorner_nodes + tot_edges_internal_nodes + tot_faces_internal_nodes + 1;
	
	for(int k=1; k<ngl-1; k++) {
	    zeta = lgl.ksi[k];
	    
	    for(int j=1; j<ngl-1; j++) {
		eta = lgl.ksi[j];
		
		for(int i=1; i<ngl-1; i++) {
		    xi = lgl.ksi[i];
		
		    COORDS_HO[ip][0] = xa*(1 - xi)*(1 - eta)*(1 - zeta)*0.125 +	\
			xb*(1 + xi)*(1 - eta)*(1 - zeta)*0.125 +	\
			xc*(1 + xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			xd*(1 - xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			xe*(1 - xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			xf*(1 + xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			xg*(1 + xi)*(1 + eta)*(1 + zeta)*0.125 +	\
			xh*(1 - xi)*(1 + eta)*(1 + zeta)*0.125;
		    
		    COORDS_HO[ip][1]  = ya*(1 - xi)*(1 - eta)*(1 - zeta)*0.125 + \
			yb*(1 + xi)*(1 - eta)*(1 - zeta)*0.125 +	\
			yc*(1 + xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			yd*(1 - xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			ye*(1 - xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			yf*(1 + xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			yg*(1 + xi)*(1 + eta)*(1 + zeta)*0.125 +	\
			yh*(1 - xi)*(1 + eta)*(1 + zeta)*0.125;
		    
		    COORDS_HO[ip][2] = za*(1 - xi)*(1 - eta)*(1 - zeta)*0.125 + \
			zb*(1 + xi)*(1 - eta)*(1 - zeta)*0.125 +	\
			zc*(1 + xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			zd*(1 - xi)*(1 + eta)*(1 - zeta)*0.125 +	\
			ze*(1 - xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			zf*(1 + xi)*(1 - eta)*(1 + zeta)*0.125 +	\
			zg*(1 + xi)*(1 + eta)*(1 + zeta)*0.125 +	\
			zh*(1 - xi)*(1 + eta)*(1 + zeta)*0.125;
		    
		    fprintf(fileidHO_vol, " %d %f %f %f\n", ip, COORDS_HO[ip][0], COORDS_HO[ip][1], COORDS_HO[ip][2]);
				
		    //IPrenumbered = IP;
		    //CONN_renumbered(IEL, iconn) = IP;
				
		    ip = ip + 1;
		    //iconn = iconn + 1;
		}
	    }
	}
    }
    
    
    return 0;
    if( ip != nnodes_linear+tot_edges_internal_nodes+tot_faces_internal_nodes+1) {
	printf(" !!!! ERROR in BUILD_CONN.c at line ~857 \n !!!! ip %d != %d nnodes_linear+tot_edges_internal_nodes+tot_faces_internal_nodes+1\n", ip,  nnodes_linear+tot_edges_internal_nodes+tot_faces_internal_nodes+1);
	printf(" !!!! The program will EXIT now\n");
	exit(1);
    }
  
    fclose(fileid);
    fclose(fileidHO_edges);
    fclose(fileidHO_faces);
    fclose(fileidHO_vol);
    
    //int NPcurrent = ip;
    
    printf(" # POPULATE GRID with SPECTRAL NODES............................ DONE\n");
    printf(" #------------------------------------------------------------------#\n");
    
    return 0;
    
}
