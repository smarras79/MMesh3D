/*
 * This function builds the connectivity matrix:
 *
 * Simone Marras, Aug 2012
 */

#include <stdlib.h>

#include "BUILD_CONN.h"
#include "GRID2CONN.h"
#include "GLOBAL_VARS.h"
#include "MEMORY.h"
#include "MYDEFINE.h"
#include "NRUTIL.h"
#include "SURFACES.h"

//BDY CODES for boundary edges:
#define TOP_FLG   111;
#define BOTT_FLG  222;
#define LEFT_FLG  222;
#define RIGHT_FLG 222;

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
    int iedge;
    int NORTH, SOUTH, EAST, WEST, BOTT, TOP;
    
    MEMORY_ALLOCATE(6);
    
    for (iel = 0; iel<nelem; iel++)
	{
	/*
	 * Edges bottom face:
	 */
	iedge = 0;
	conn_edge_el[iel][iedge][0] = CONN[iel][0];
	conn_edge_el[iel][iedge][1] = CONN[iel][1];
	iedge = 1;
	conn_edge_el[iel][iedge][0] = CONN[iel][1];
	conn_edge_el[iel][iedge][1] = CONN[iel][2];
	iedge = 2;
	conn_edge_el[iel][iedge][0] = CONN[iel][2];
	conn_edge_el[iel][iedge][1] = CONN[iel][3];
	iedge = 3;
	conn_edge_el[iel][iedge][0] = CONN[iel][3];
	conn_edge_el[iel][iedge][1] = CONN[iel][0];
	//Edges top face
	iedge = 4;
	conn_edge_el[iel][iedge][0] = CONN[iel][4];
	conn_edge_el[iel][iedge][1] = CONN[iel][5];
	iedge = 5;
	conn_edge_el[iel][iedge][0] = CONN[iel][5];
	conn_edge_el[iel][iedge][1] = CONN[iel][6];
	iedge = 6;
	conn_edge_el[iel][iedge][0] = CONN[iel][6];
	conn_edge_el[iel][iedge][1] = CONN[iel][7];
	iedge = 7;
	conn_edge_el[iel][iedge][0] = CONN[iel][7];
	conn_edge_el[iel][iedge][1] = CONN[iel][4];
	
	//Vertical edges
	iedge = 8;
	conn_edge_el[iel][iedge][0] = CONN[iel][0];
	conn_edge_el[iel][iedge][1] = CONN[iel][4];
	iedge = 9;
	conn_edge_el[iel][iedge][0] = CONN[iel][1];
	conn_edge_el[iel][iedge][1] = CONN[iel][5];
	iedge = 10;
	conn_edge_el[iel][iedge][0] = CONN[iel][2];
	conn_edge_el[iel][iedge][1] = CONN[iel][6];
	iedge = 11;
	conn_edge_el[iel][iedge][0] = CONN[iel][3];
	conn_edge_el[iel][iedge][1] = CONN[iel][7];
	
	/*
	 * Local faces node connectivity:
	 * i.e. what nodes belong to a given local face in iel:
	 */
	SOUTH = 0; //south
	conn_face_el[iel][SOUTH][0] = CONN[iel][0];
	conn_face_el[iel][SOUTH][1] = CONN[iel][1];
	conn_face_el[iel][SOUTH][2] = CONN[iel][5];
	conn_face_el[iel][SOUTH][3] = CONN[iel][4];
	NORTH = 2; //north
	conn_face_el[iel][NORTH][0] = CONN[iel][2];
	conn_face_el[iel][NORTH][1] = CONN[iel][3];
	conn_face_el[iel][NORTH][2] = CONN[iel][7];
	conn_face_el[iel][NORTH][3] = CONN[iel][6];
	EAST = 1; //east
	conn_face_el[iel][EAST][0] = CONN[iel][1];
	conn_face_el[iel][EAST][1] = CONN[iel][2];
	conn_face_el[iel][EAST][2] = CONN[iel][6];
	conn_face_el[iel][EAST][3] = CONN[iel][5];
	WEST = 3; //west
	conn_face_el[iel][WEST][0] = CONN[iel][3];
	conn_face_el[iel][WEST][1] = CONN[iel][0];
	conn_face_el[iel][WEST][2] = CONN[iel][4];
	conn_face_el[iel][WEST][3] = CONN[iel][7];
	BOTT = 4; //bottom
	conn_face_el[iel][BOTT][0] = CONN[iel][3];
	conn_face_el[iel][BOTT][1] = CONN[iel][2];
	conn_face_el[iel][BOTT][2] = CONN[iel][1];
	conn_face_el[iel][BOTT][3] = CONN[iel][0];
	TOP = 5; //top
	conn_face_el[iel][TOP][0] = CONN[iel][4];
	conn_face_el[iel][TOP][1] = CONN[iel][5];
	conn_face_el[iel][TOP][2] = CONN[iel][6];
	conn_face_el[iel][TOP][3] = CONN[iel][7];
	
	/*for (int iface=0; iface<6; iface++)
	    {
		iedge = 0;
		conn_face_edge_el[iel][iface][iedge][0] = conn_face_el[iel][iface][1];
		conn_face_edge_el[iel][iface][iedge][1] = conn_face_el[iel][iface][2];
		conn_face_edge_el[iel][iface][iedge][2] = iel;
		
		iedge = 1;
		conn_face_edge_el[iel][iface][iedge][0] = conn_face_el[iel][iface][2];
		conn_face_edge_el[iel][iface][iedge][1] = conn_face_el[iel][iface][3];
		conn_face_edge_el[iel][iface][iedge][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge][iface][iel][conn_face_edge_el[iel][iface][iedge][0]][conn_face_edge_el[iel][iface][iedge][1]];
		
		iedge = 2;
		conn_face_edge_el[iel][iface][iedge][0] = conn_face_el[iel][iface][3];
		conn_face_edge_el[iel][iface][iedge][1] = conn_face_el[iel][iface][4];
		conn_face_edge_el[iel][iface][iedge][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge][iface][iel][conn_face_edge_el[iel][iface][iedge][0]][conn_face_edge_el[iel][iface][iedge][1]];
		
		iedge = 3;
		conn_face_edge_el[iel][iface][iedge][0] = conn_face_el[iel][iface][4];
		conn_face_edge_el[iel][iface][iedge][1] = conn_face_el[iel][iface][1];
		conn_face_edge_el[iel][iface][iedge][2] = iel;
		//fprintf[' edge %d of face %d of element %d has points [%d][%d]\n'][iedge][iface][iel][conn_face_edge_el[iel][iface][iedge][0]][conn_face_edge_el[iel][iface][iedge][1]];
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
    int conn_face_el_sort[nelem][6][4];
    int AUXI[5];
    for (int i=0; i<nelem; i++) {
	for (int j=0; j<6; j++) {
	    for (int k=0; k<4; k++) {
		int kk = k + 1;
		AUXI[kk] = conn_face_el[i][j][k];
	    }
	    isort(4, AUXI);
	    for (int k=0; k<4; k++) {
		int kk = k + 1;
		conn_face_el_sort[i][j][k] = AUXI[kk];
		    printf(" %d\n", conn_face_el_sort[i][j][k]);
	    } printf("\n");
	}
    }
    //END sorting 
    

    iface = 0;
    for(iel = 0; iel<nelem; iel++) {
	for(int jel = 0; jel<nelem; jel++) {
		

	    if(     conn_face_el_sort[iel][BOTT][0] == conn_face_el_sort[jel][TOP][0] && \
		    conn_face_el_sort[iel][BOTT][1] == conn_face_el_sort[jel][TOP][1] && \
		    conn_face_el_sort[iel][BOTT][2] == conn_face_el_sort[jel][TOP][2] && \
		    conn_face_el_sort[iel][BOTT][3] == conn_face_el_sort[jel][TOP][3])
		{			
		    FACE_in_ELEM[iel][BOTT][0] = iel;
		    FACE_in_ELEM[iel][BOTT][1] = jel;
				
		    FACE_in_ELEM[jel][TOP][0] = jel;
		    FACE_in_ELEM[jel][TOP][1] = iel;
				
		    //fprintf(' SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n'][BOTT][iel][TOP][jel];
				
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
			
		    //fprintf(' SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n'][EAST][iel][WEST][jel];
			
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
			
		    //fprintf(' SHARED FACE:  face %d of ELEMENT %d is shared with face %d of ELEMENT %d\n', SOUTH, iel, NORTH, jel];
		    iface = iface + 1;
		}
	}
    }
		    
    int nint_faces        = iface - 1;
    int nfaces            = nint_faces + nbdy_faces;
    int face_direction[6] = {SOUTH, EAST, NORTH, WEST, BOTT, TOP};

    /*----------------------------------------------------------------------------------------------------
     * Store connectivity of internal faces: 
     * start from the boundary faces and append the internal faces afterwards:
     *----------------------------------------------------------------------------------------------------*/
    //1) Store bdy faces in CONN_FACE(1:BDY_FACES, 1:4)
    int CONN_FACE[nfaces][4];
    int CONN_FACE_all[nfaces][4];
    
    for (int i=0; i<nbdy_faces; i++) {
	for (int j=0; j<4; j++) {
	    CONN_FACE[i][j] = CONN_BDY_FACES[i][j];
	}
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
    int iel;
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
    
    int nface_all = iface_all - 1;
    int REPEATED_auxi[nfaces][2];
    int REPEATED_index[nface_all];
    int FACE_MULTIPLICITY_auxi[nface_all];
    for (int i=0; i<nfaces; i++) {
	REPEATED_auxi[i][0] = -1;
	REPEATED_auxi[i][1] = -1;
    }    
    for (int i=0; i<nface_all; i++) {
	FACE_MULTIPLICITY_auxi[i] = 0;
    }

    int iface_repeated = 0;
    int ifac           = 1;
    int krepeated      = 0;
    for (int i=0; i<nface_all; i++) {
	int multiplicity = 0;
	for (int j=i; j<nface_all; j++) {
		
	    if (j != i) {
		if ((CONN_FACE_all[i][0] == CONN_FACE_all[j][0] &&	\
		     CONN_FACE_all[i][1] == CONN_FACE_all[j][1] &&	\
		     CONN_FACE_all[i][2] == CONN_FACE_all[j][2] &&	\
		     CONN_FACE_all[i][3] == CONN_FACE_all[j][3]) ) {
		    
		    multiplicity = multiplicity  + 1;
		    krepeated = krepeated + 1;
		    REPEATED_index[krepeated] = j;
		    FACE_MULTIPLICITY_auxi[i] = multiplicity + 1;
		    
		    //printf('k=%d, REPEATED_index %d -> (%d %d %d %d) repeated %d times\n', krepeated, REPEATED_index(krepeated), CONN_FACE_tmp(j, 1), CONN_FACE_tmp(j, 2), CONN_FACE_tmp(j, 3), CONN_FACE_tmp(j, 4), multiplicity);
		}
	    }
	}
    }

    
    /*
     * Set repeated entries of CONN_EDGE_tmp to -1:
     */
    int nrepeated = krepeated;
    for (int i=0; i<nrepeated; i++) {
	int irepeated_index = REPEATED_index[i];
	
	CONN_FACE_all[irepeated_index][0] = -1;
	CONN_FACE_all[irepeated_index][1] = -1;
	CONN_FACE_all[irepeated_index][2] = -1;
	CONN_FACE_all[irepeated_index][3] = -1;
    }

    
    /*
     * STORE EACH FACE into CONN_FACE(1:`, 1:4):
     */
    ifac = 0;
    int FACE_MULTIPLICITY[nface_all]; //overallocated
    for (int iface_all=0; iface_all<nface_all; iface_all++) {
	if (CONN_FACE_all[iface_all][0] > 0) {
	    
	    CONN_FACE[ifac][0] = CONN_FACE_all[iface_all][0];
	    CONN_FACE[ifac][1] = CONN_FACE_all[iface_all][1];
	    CONN_FACE[ifac][2] = CONN_FACE_all[iface_all][2];
	    CONN_FACE[ifac][3] = CONN_FACE_all[iface_all][3];
	    /*
	      FACE_MULTIPLICITY[ifac] = FACE_MULTIPLICITY_auxi[iface_all];
	    */
	    // FACE_in_ELEM[ifac][0] = FACE_in_ELEM[iface_all][0];
	    // FACE_in_ELEM[ifac][1] = FACE_in_ELEM[iface_all][1];
	    
	    ifac = ifac + 1;
	}
    }
    
    printf(" NELEM %d. nface_all %d %d %d", nelem, nface_all, ifac, nfaces);
    exit(1);
    /*for (int iface=0; iface<nfaces; iface++) {
	printf(" CONN_FACE(%d,1:4) = %d %d %d %d, has multiplicity %d\n", iface, CONN_FACE[iface][0], CONN_FACE[iface][1], CONN_FACE[iface][2], CONN_FACE[iface][3], FACE_MULTIPLICITY[iface]);
	}*/

    /*--------------------------------------------------------------------------
     * Populate FACE_LtoG(1:NEL,1:6) OK
     *--------------------------------------------------------------------------*
    for (int iface=0; iface<nfaces; iface++) {
	for (int iel=0; iel<nelem; iel++) {
	    for (int iface=0; iface<6; iface++) {
		
		if ( (  CONN_FACE[iface][0] == conn_face_el[iel][iface][0] && \
			CONN_FACE[iface][1] == conn_face_el[iel][iface][1] && \
			CONN_FACE[iface][2] == conn_face_el[iel][iface][2] && \
			CONN_FACE[iface][3] == conn_face_el[iel][iface][3]) ) {
			    
			    FACE_LtoG[iel][iface] = iface;
			    //%fprintf('  --- FACE_LtoG(%d,%d) = %d -> (%d %d %d %d) \n', iel, iface, FACE_LtoG(iel, iface), conn_face_el_sort[iel][iface][0], \
			    //	%conn_face_el_sort[iel][iface][1], \
			    //	%conn_face_el_sort[iel][iface][2], \
			    //	%conn_face_el_sort[iel][iface][3]);
		}
	    }
	}
	}*/
    
    
    return 0;
}

int ADD_HIGH_ORDER_NODES(void)
{
    /*--------------------------------------------------------------------------
     * Auxiliary elemental parameters/counters:
     *--------------------------------------------------------------------------*
    int ncorner_nodes            = 8;
    int nedge_internal_nodes     = (ngl-2);
    int nface_internal_nodes     = (ngl-2)*(ngl-2);
    int nvolume_internal_nodes   = (ngl-2)*(ngl-2)*(ngl-2);
    int tot_edges_internal_nodes = nedge_internal_nodes*12;
    int tot_faces_internal_nodes = nface_internal_nodes*6;

    int iedge = 1;
    int ip = nnodes + 1; //we start populating from the low order numbering
    
    printf(" #------------------------------------------------------------------#\n");
    printf(" # POPULATE GRID with SPECTRAL NODES............................\n");

    printf(" nnodes = %d\n", nnodes);

    /*--------------------------------------------------------------------------
     * Build high order grid points on every edges
     --------------------------------------------------------------------------*
    for(int iel = 0; iel<nelem; iel++)
	{
	    int iconn = 9;
	    for(int iedg = 0; iedg<12; iedg++)
		{
		    iedge = EDGE_LtoG[iel][iedg];
		    //printf(' iedge %d of iel %d is GLOBA iedge %d\n'][iedg, iel, iedge);
		    if (EDGE_FLAG[iedge] < 0) {
			EDGE_FLAG[iedge] = 1;
			
			int ip1 = min(CONN_EDGE[iedge][1], CONN_EDGE[iedge][2]);
			int ip2 = max(CONN_EDGE[iedge][1], CONN_EDGE[iedge][2]);
			
			double x1  = COORDS[ip1][0], y1  = COORDS[ip1][1], z1  = COORDS[ip1][2];
			double x2  = COORDS[ip2][0], y2  = COORDS[ip2][1], z2  = COORDS[ip2][2];
			
			/*
			 * Plot LGL points on edges:
			 *
			for(int l = 2; l<=ngl-1; l++)
			    {
				/*
				xi = x_LGL[l];
				
				xx[l] = x1*(1.0 - xi)*0.5 + (1.0 + xi)*x2*0.5;
				yy[l] = y1*(1.0 - xi)*0.5 + (1.0 + xi)*y2*0.5;
				zz[l] = z1*(1.0 - xi)*0.5 + (1.0 + xi)*z2*0.5;
				COORDS[IP][2) = xx[l];
				COORDS[IP][3) = yy[l];
				COORDS[IP][4) = zz[l];
				EDGE_POINT_CONN[iedge][l] = IP;
				EDGE_CONN[iedge][l] = IP;
				
				IP = IP + 1; //Initialized to highest low order value of npoin.
				iconn = iconn + 1;
				*
			    }
		    }//	end if
		} //end iedg
	} //end iel
    int NPcurrent = ip;
    

    printf(" # POPULATE GRID with SPECTRAL NODES............................ DONE\n");
    printf(" #------------------------------------------------------------------#\n");
    */
    return 0;
    
}
