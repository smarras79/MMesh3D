#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Variables declarations:
#include "MYDEFINE.h"

/*
 * Scalars
 */
double	 						\
xlength,						\
    ylength,						\
    zlength,						\
    xmin,						\
    xmax,						\
    ymin,						\
    ymax,						\
    zmin,					       	\
    zmax,				       		\
    xmin_l,			       			\
    xmax_l,						\
    ymin_l,						\
    ymax_l,						\
    zmin_l,						\
    zmax_l,						\
    xmin_obs,						\
    xmax_obs,						\
    ymin_obs,						\
    ymax_obs,						\
    zmin_obs,						\
    zmax_obs,						\
    xc,							\
    yc,							\
    zc,							\
    dx,							\
    dy,							\
    dz,							\
    known_x,						\
    known_y,						\
    known_lat,						\
    known_lon,						\
    deltaLon,						\
    deltaLat,						\
    scalefactor;

int									\
nnodesx_obs,							\
    nnodesy_obs,							\
    nnodesz_obs,							\
    wW,									\
    wE,									\
    wN,									\
    wS,									\
    obs_wW,								\
    obs_wE,								\
    obs_wN,								\
    obs_wS,								\
    nnodes_g,								\
    nnodes,								\
    nnodes_sem,								\
    nnodesx,								\
    nnodesy,								\
    nnodesz,								\
    nx,									\
    ny,									\
    nz,									\
    nelem,								\
    nelem_g,								\
    nfaces,								\
    nedges,								\
    nbdy_faces,								\
    nbdy_edges,								\
    nelx,								\
    nely,								\
    nelz,								\
    ncol,								\
    ncol_g,								\
    narray,								\
    ngl,								\
    nglx,								\
    ngly,								\
    nglz,								\
    nop,								\
    nopx,								\
    nopy,								\
    nopz,								\
    nboun,								\
    nboun_g,								\
    elem,								\
    EL_NODES,								\
    i,									\
    j,									\
    k,									\
    iel,								\
    ipoin,								\
    status,								\
    readgeo_flg,							\
    line_cntr,								\
    cntr,								\
    count,								\
    len,								\
    isigned,								\
    endian,								\
    mpiprocs,								\
    wordsize,								\
    lMETIS,								\
    lCART,								\
    lread_external_grid;

/*
 * Pointers
 */
double									\
*rarray,								\
    legendre,								\
    dlegendre,								\
    *lagrange,								\
    *dlagrange,								\
    **COORDS,								\
    *COORDS1d,								\
    **BDY_COORDS,							\
    **BDYFLAG,								\
    **PHI,								\
    *parameters,							\
    *xgl,								\
    *wgl;

double						\
***x,						\
    ***y,					\
    ***z;

int								\
**CONN,								\
    **CONN_FACE,						\
    **CONN_EDGE,						\
    **CONN_BDY_EDGES,						\
    **CONN_BDY_FACES,						\
    ***conn_edge_el,						\
    ***conn_face_el,						\
    ***conn_face_edge_el,					\
    ***FACE_in_ELEM,						\
    **FACE_LtoG,						\
    **EDGE_LtoG,						\
    *ELTYPE;

char								    	\
*inputfile,								\
    *external_grid_file_name,						\
    *print_alya,							\
    *print_vtk,								\
    *print_gmsh,							\
    *fname,								\
    *input_inp,								\
    *outfile_msh_vtk,							\
    *outfile_msh_gmsh,							\
    *outfile_msh_alya,							\
    *problem[PROB_ENTRIES];

/*
 * Static arrays (used for small arrays only) 
 */
char									\
nodes  [6],								\
    nodesx [6],								\
    nodesy [6],								\
    nodesz [6],								\
    elorder[3],								\
    nel[5],								\
//NOTE: if you want to store 3 digits, you need 4 spaces because C stores "/0" as last. If you don't, you will get the Abort Trap error! 
    mpiprocess[5];

int						\
INPUTVariables[MAX_INPUTS],			\
    grid_type[3];

/* Notes on grid_type:
 * Vector of 3 elements: each element can be either any integer that defines a type of
 * grid points distribution. This should be defined by the user in the input.inp file.
 *	1 <-- uniform grid with equispaced grid points
 *	2 <-- exponential grid points distribution from left to right [x(i)=exp(i)-1]
 *	3 <-- like 2, but from right to left (TO BE CODED YET)
 *
 *	Ex.:
 *	grd_type[2, 1, 1] means: that in x we want an exponential distribution,
 *				 while in y and z a uniform grid.
 * Lattice coordinates for the INTERPOLATE function (PHI[1:m][1:n], where m,n are user inputs)
 */

/*
 * FILE pointers
 */
FILE						\
*input_id,					\
    *vtk_file_id,				\
    *gmsh_file_id;

#endif
