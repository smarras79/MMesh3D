/*
   p4est is a C library to manage a collection (a forest) of multiple
   connected adaptive quadtrees or octrees in parallel.

   Copyright (C) 2010 The University of Texas System
   Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

   p4est is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   p4est is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with p4est; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// #define VTK_OUTPUT

/*#ifndef P4_TO_P8
#define P4_TO_P8
#endif
*/
#include <sc.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef P4_TO_P8
#include <p8est_bits.h>
#include <p8est_connectivity.h>
//#include <p8est_connrefine.h> //Sohail Reddy (08/05/2020): comment out to use with new p4est
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>
#else
#include <p4est_bits.h>
#include <p4est_connectivity.h>
//#include <p4est_connrefine.h> //Sohail Reddy (08/05/2020): comment out to use with new p4est
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#endif


#include "P4EST_API.h"

// define single/double precision
#ifdef SINGLE
#define real float
#else
#define real double
#endif

#ifdef P4_TO_P8
#define FUN_PREFIX p8esttonuma_
#else
#define FUN_PREFIX p4esttonuma_
#endif

#if ibmcompiler
#define FUN_POSTFIX
#else
#define FUN_POSTFIX _
#endif

#define P2N_APPEND_3(a, b, c) a##b##c
#define P2N_APPEND_EXPAND_3(a, b, c) P2N_APPEND_3(a, b, c)
#define P2N_FUN(X) P2N_APPEND_EXPAND_3(FUN_PREFIX, X, FUN_POSTFIX)

typedef struct
{
  int iref;
  int oref;
} user_data_t;

typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpisize;
  int mpirank;
} mpi_context_t;

/* Defines a structure which exports data to NUMA*/
typedef struct
{
  int nop;
  p4est_locidx_t npoin;
  p4est_locidx_t npoin_dg;
  p4est_locidx_t nface;
  p4est_locidx_t nface_dg;
  p4est_locidx_t nelem;

  p4est_locidx_t *vc_el_type; // list of sponge types

  int num_nbh;
  real *coord_dg;
  real *coord_cg;

  p4est_locidx_t *intma_table;
  int *NC_face;
  int *NC_edge;
  int num_send_recv_total;
  int *num_send_recv;
  int *nbh_proc;
  int *nbh_send_recv;
  int *nbh_send_recv_multi;
  int *nbh_send_recv_half;

  p4est_locidx_t *EToNC;
  p4est_locidx_t nNC;

  p4est_locidx_t nbsido;
  p4est_locidx_t nboun;
  p4est_locidx_t *bsido;
  p4est_locidx_t *face;
  p4est_locidx_t *face_type;
  p4est_locidx_t *level_list;
  p4est_locidx_t *plist;
} p4esttonuma_t;

/* Defines a structure to hold BC information from a file*/
typedef struct
{
  p4est_topidx_t num_bc;
  p4est_topidx_t *bc_to_vertex;
  p4est_topidx_t *bc_physical_type;
  p4est_topidx_t *bc_el_label;
  p4est_topidx_t *vc_el_type;
  p4est_topidx_t *vc_el_label;
  p4est_locidx_t nel_vc_sponge; // numbers of elements with sponge
} p4est_bc_info_t;

/* ------------------------------------- */
/* Save previous forest and connectivity */
/* ------------------------------------- */
static p4est_t *stored_p4est = NULL;
static int refine_level = 0;
static p4est_bc_info_t *bc_data;

#ifdef P4_TO_P8
extern mpi_context_t mpi_context;
#else
mpi_context_t mpi_context;
#endif

/**
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 */
static void fill_coordinates(int Nrp, real *r, p4est_t *p4est,
                             p4est_lnodes_t *lnodes, p4esttonuma_t *p2n)
{
  int i, c; /* We use plain int for small loops. */
  int m, l;
  double vxyz[P4EST_CHILDREN][3]; /* We embed the 2D vertices into 3D space. */
  p4est_topidx_t tt;              /* Connectivity variables have this type. */
  p4est_locidx_t k, q, Q;         /* Process-local counters have this type. */
  p4est_tree_t *tree;             /* Pointer to one octree */
  p4est_quadrant_t *quad, node;
  sc_array_t *tquadrants; /* Quadrant array for one tree */

  /* Allocate coordinates in DG storage */
  p2n->coord_dg = (real *)malloc(sizeof(real) * 3 * p2n->npoin_dg);
  p2n->coord_cg = (real *)malloc(sizeof(real) * 3 * p2n->npoin);

  /* Go over all elements and compute coordinated of nodal points*/
  for (tt = p4est->first_local_tree, k = 0; tt <= p4est->last_local_tree; ++tt)
  {

    tree = p4est_tree_array_index(p4est->trees, tt); /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t)tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k)
    {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index(tquadrants, q);

      for (i = 0; i < P4EST_CHILDREN; ++i)
      {
        p4est_quadrant_corner_node(quad, i, &node);

/* Transform per-tree reference coordinates into physical space. */
#ifdef P4_TO_P8
        p4est_qcoord_to_vertex(p4est->connectivity, tt, node.x, node.y, node.z,
                               vxyz[i]);
#else
        p4est_qcoord_to_vertex(p4est->connectivity, tt, node.x, node.y,
                               vxyz[i]);
#endif
      }

/* Loop over element degrees of freedom to fill the coordinates */
/* coord(3,1:npoin): x,y,z-coordinates of each point */
#ifdef P4_TO_P8
      for (int n = 0; n < Nrp; ++n)
#endif
        for (m = 0; m < Nrp; ++m)
          for (l = 0; l < Nrp; ++l)
          {
            real w[P4EST_CHILDREN];
#ifdef P4_TO_P8
            p4est_locidx_t offset = n * Nrp * Nrp + m * Nrp + l;
            w[0] = (1 - r[l]) * (1 - r[m]) * (1 - r[n]);
            w[1] = (1 + r[l]) * (1 - r[m]) * (1 - r[n]);
            w[2] = (1 - r[l]) * (1 + r[m]) * (1 - r[n]);
            w[3] = (1 + r[l]) * (1 + r[m]) * (1 - r[n]);
            w[4] = (1 - r[l]) * (1 - r[m]) * (1 + r[n]);
            w[5] = (1 + r[l]) * (1 - r[m]) * (1 + r[n]);
            w[6] = (1 - r[l]) * (1 + r[m]) * (1 + r[n]);
            w[7] = (1 + r[l]) * (1 + r[m]) * (1 + r[n]);
            const p4est_locidx_t dg_index = Nrp * Nrp * Nrp * k + offset;
#else
          p4est_locidx_t offset = m * Nrp + l;
          w[0] = (1 - r[l]) * (1 - r[m]);
          w[1] = (1 + r[l]) * (1 - r[m]);
          w[2] = (1 - r[l]) * (1 + r[m]);
          w[3] = (1 + r[l]) * (1 + r[m]);
          const p4est_locidx_t dg_index = Nrp * Nrp * k + offset;
#endif
            const p4est_locidx_t cg_index = lnodes->element_nodes[dg_index];
            real temp;
            for (i = 0; i < 3; ++i)
            {
              temp = 0;
              for (c = 0; c < P4EST_CHILDREN; c++)
                temp += w[c] * vxyz[c][i];
              temp /= P4EST_CHILDREN;
              p2n->coord_dg[3 * dg_index + i] = temp;
              p2n->coord_cg[3 * cg_index + i] = temp;
            }
          }
    }
  }
}
/* A function to initialize user data for each quadrant- here used to
   initialize refinement flag to 0*/
static void init_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  user_data_t *data = (user_data_t *)quadrant->p.user_data;

  data->iref = 0;
}

/* Function returns a refine flag for all elements below maximum refinement
   level
   enforcing uniform domain refinement*/
static int refine_uniform_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *quadrant)
{
  return (int)quadrant->level < refine_level;
}

/*Function returns a refinement flag which is stored in quadrant user_data*/
static int refine_list_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t *quadrant)
{
  user_data_t *data = (user_data_t *)quadrant->p.user_data;

  return (int)(data->iref > 0);
}

/*Function returns a coarsening flag which is stored in quadrant user_data*/
static int coarsen_list_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *q[])
{
  user_data_t *data;

  p4est_locidx_t i;
  for (i = 0; i < P4EST_CHILDREN; i++)
  {
    data = (user_data_t *)q[i]->p.user_data;
    if (data->iref >= 0)
      return 0;
  }

  return 1;
}
/* Cubed shell */
#ifndef P4_TO_P8

static p4est_connectivity_t *p4est_connectivity_new_cubed_sphere(void)
{
  const p4est_topidx_t num_vertices = 8;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_ctt = 0;
  const double vertices[8 * 3] = {
      -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1,
      -1, -1, 1,  1, -1, 1,  -1, 1, 1,  1, 1, 1,
  };
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
      0, 2, 1, 3, 2, 6, 3, 7, 0, 4, 2, 6, 4, 5, 6, 7, 0, 1, 4, 5, 1, 3, 5, 7,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
      4, 1, 2, 5, 0, 3, 2, 5, 0, 3, 4, 1, 2, 5, 4, 1, 2, 5, 0, 3, 4, 1, 0, 3,
  };
  const int8_t tree_to_face[6 * 4] = {
      2, 0, 0, 2, 1, 3, 3, 1, 2, 0, 0, 2, 1, 3, 3, 1, 2, 0, 0, 2, 1, 3, 3, 1,
  };

  return p4est_connectivity_new_copy(num_vertices, num_trees, 0, vertices,
                                     tree_to_vertex, tree_to_tree, tree_to_face,
                                     NULL, &num_ctt, NULL, NULL);
}

#else

static p4est_connectivity_t *p8est_connectivity_new_cubed_sphere(void)
{
  const p4est_topidx_t num_vertices = 16;
  const p4est_topidx_t num_trees = 6;
  const double vertices[16 * 3] = {
      +1, -1, -1, // 0
      +1, +1, -1, // 1
      +1, -1, +1, // 2
      +1, +1, +1, // 3
      -1, -1, +1, // 4
      -1, +1, +1, // 5
      -1, -1, -1, // 6
      -1, +1, -1, // 7
      +2, -2, -2, // 8
      +2, +2, -2, // 9
      +2, -2, +2, // 10
      +2, +2, +2, // 11
      -2, -2, +2, // 12
      -2, +2, +2, // 13
      -2, -2, -2, // 14
      -2, +2, -2, // 15
  };
  //       0---1
  //       | 3 |
  //   0---6---7---1
  //   | 5 | 2 | 4 |
  //   2---4---5---3
  //       | 1 |
  //       2---3
  //       | 0 |
  //       0---1
  const p4est_topidx_t tree_to_vertex[6 * 8] = {
      0, 1, 2, 3, 8,  9,  10, 11, // 0
      2, 3, 4, 5, 10, 11, 12, 13, // 1
      4, 5, 6, 7, 12, 13, 14, 15, // 2
      6, 7, 0, 1, 14, 15, 8,  9,  // 3
      5, 3, 7, 1, 13, 11, 15, 9,  // 4
      2, 4, 0, 6, 10, 12, 8,  14, // 5
  };

  p4est_connectivity_t *conn =
      p4est_connectivity_new(num_vertices, num_trees, 0, 0, 0, 0);
  for (int tree = 0; tree < conn->num_trees; ++tree)
  {
    // Fill tree_to_vertex
    for (int c = 0; c < P4EST_CHILDREN; ++c)
      conn->tree_to_vertex[P4EST_CHILDREN * tree + c] =
          tree_to_vertex[P4EST_CHILDREN * tree + c];

    // Fill with data to make valid connectivity
    for (int face = 0; face < P4EST_FACES; ++face)
    {
      conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
      conn->tree_to_face[P4EST_FACES * tree + face] = face;
    }
  }

  // Fill vertices
  for (int n = 0; n < 3 * num_vertices; ++n)
    conn->vertices[n] = vertices[n];

  // Compute real tree_to_* fields and complete (edge and) corner fields.
  p4est_connectivity_complete(conn);

  return conn;
}

#endif

/* Start p4est */
#ifdef P4_TO_P8

#ifdef ibmcompiler
void p4est_start
#else
void p4est_start_
#endif
    (int *p4est_log_level)
{
  int mpiret;
  mpi_context_t *mpi = &mpi_context;
  /* initialize MPI and p4est internals */
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size(mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI(mpiret);
  mpiret = sc_MPI_Comm_rank(mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI(mpiret);

  sc_init(mpi->mpicomm, 0, 1, NULL, *p4est_log_level);
  p4est_init(NULL, *p4est_log_level);
}

#ifdef ibmcompiler
void p4est_stop
#else
void p4est_stop_
#endif
    ()
{
  p4est_t *p4est = stored_p4est;
  if (p4est) // don't know why this could be null !?
  {
    p4est_connectivity_t *connectivity = p4est->connectivity;
    p4est_destroy(p4est);
    if (connectivity)
      p4est_connectivity_destroy(connectivity);
    connectivity = NULL;
  }
  p4est = NULL;
  stored_p4est = NULL;
  //   sc_finalize(); //causes crash un-deallocated memory ?
}

#endif

static int FACE_LEN;
static int ORIENT;

static void quadrant_init_replace(p4est_t *p4est, p4est_topidx_t which_tree,
                                  int num_outgoing,
                                  p4est_quadrant_t *outgoing[],
                                  int num_incoming,
                                  p4est_quadrant_t *incoming[])
{
  SC_ASSERT(num_outgoing == 1);
  SC_ASSERT(num_incoming == P4EST_CHILDREN);

  const user_data_t *outdata = (user_data_t *)outgoing[0]->p.user_data;

  int c = 0;
  for (c = 0; c < P4EST_CHILDREN; ++c)
  {
    user_data_t *indata = (user_data_t *)incoming[c]->p.user_data;
    indata->iref = outdata->iref - 1;
  }
}

// Function compares two integers
static int compare(const void *a, const void *b)
{
  return (*(int *)a - *(int *)b);
}

// Function checks which bc is prescribed in the line of text
static int which_bc(char *line)
{

  const int num_bc_types = 13;
  const char *bc_list_label[13] = {"FREE_SLIP",   "NO_SLIP",  "T_DIRICHLET",
                                  "S_DIRICHLET", "OUTFLOW",  "U_SPONGE",
				   "T_SPONGE",    "S_SPONGE", "TS_ICE", "PERIODIC", "TOP_SPONGE", "INFLOW", "OUTFLOW"};
  const int bc_list_condition[13] = {4, 5, 50, 500, 6, 1, 10, 100, 770, 3, 6, 5, 6};

  int bc = 0;

  for (int i = 0; i < num_bc_types; i++)
  {
    if (strstr(line, bc_list_label[i]))
    {
      bc += bc_list_condition[i];
      // return bc_list_condition[i];
    }
  }
  return bc;
}

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
static char *p4est_connectivity_getline_upper(FILE *stream)
{
  char *line = P4EST_ALLOC(char, 1024), *linep = line;
  size_t lenmax = 1024, len = lenmax;
  int c;

  if (line == NULL)
    return NULL;

  for (;;)
  {
    c = fgetc(stream);
    c = toupper(c);
    if (c == EOF && linep == line)
    {
      P4EST_FREE(linep);
      return NULL;
    }

    if (--len == 0)
    {
      len = lenmax;
      lenmax *= 2;
      char *linen = P4EST_REALLOC(linep, char, lenmax);

      if (linen == NULL)
      {
        P4EST_FREE(linep);
        return NULL;
      }
      line = linen + (line - linep);
      linep = linen;
    }
    if ((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

/* This function reads the mesh file and fills bc info*/
static int p4est_bc_read_inp_stream(FILE *stream, p4est_topidx_t *num_bc,
                                    p4est_topidx_t *num_elem,
                                    p4est_topidx_t *bc_physical_type,
                                    p4est_topidx_t *bc_el_label,
                                    p4est_topidx_t *bc_to_vertex,
                                    p4est_topidx_t *vc_el_type,
                                    p4est_topidx_t *vc_el_label)
{
  int reading_elements = 0;
  int reading_elsets = 0;
  int reading_elsets_volume = 0;
  int reading_elements_volume = 0;
  int wrong_element = 0;
  int lines_read = 0, lines_free = 0;
  char *line;
  p4est_topidx_t num_physical_el = 0;
  p4est_topidx_t num_elements = 0;
  p4est_topidx_t num_elsets = 0;
  p4est_topidx_t num_elsets_volume = 0;
  p4est_topidx_t num_elements_volume = 0;
  int fill_bc_data = (bc_physical_type != NULL && bc_to_vertex != NULL);
  int boundary_condition = 0;
  int *physical_type;
  int *element_label;
  int *vc_el_label_list;
  int ibc = 0;
  int ivc = 0;

  P4EST_ASSERT((bc_physical_type == NULL && bc_to_vertex == NULL) ||
               (bc_physical_type != NULL && bc_to_vertex != NULL));

  if (fill_bc_data)
  {
    physical_type = malloc(*num_bc * sizeof(int));
    element_label = malloc(*num_bc * sizeof(int));
    vc_el_label_list = malloc(*num_elem * sizeof(int));
  }

  for (;;)
  {
    line = p4est_connectivity_getline_upper(stream);

    if (line == NULL)
    {
      break;
    }

    ++lines_read;

    /* check for control line */
    if (line[0] == '*')
    {
      reading_elements = reading_elsets = reading_elsets_volume =
          reading_elements_volume = wrong_element = 0;
      if (strstr(line, "*ELEMENT"))
      {

        if (
#ifdef P4_TO_P8
            strstr(line, "TYPE=C2D4") || strstr(line, "TYPE=CPS4") ||
            strstr(line, "TYPE=S4") ||
            strstr(line, "TYPE=CPS3") // quad elements, mostly CPS4
#else
            strstr(line, "TYPE=T3D2")
#endif
            )
        {
          if (strstr(line, "TYPE=CPS3")) // we don't want triangles
            wrong_element = 1;
          reading_elements = 1;
          ++lines_free;
          P4EST_FREE(line);
          continue;
        }
        else if (strstr(line, "TYPE=C3D8")) // just counting 3D elements
        {
          reading_elements_volume = 1;
          ++lines_free;
          P4EST_FREE(line);
          continue;
        }
      }
      else if (strstr(line, "*ELSET") &&
               strstr(line, ":BC_")) // boundary conditions (faces)
      {

        // find which bc is prescribed
        boundary_condition = which_bc(line);

        reading_elsets = 1;
        ++lines_free;
        P4EST_FREE(line);
        continue;
      }
      else if (strstr(line, "*ELSET") &&
               strstr(line, ":VC_")) // volume condition (sponge)
      {
        // find which vc is prescribed
        boundary_condition = which_bc(line);
        reading_elsets_volume = 1;
        ++lines_free;
        P4EST_FREE(line);
        continue;
      }
    }

    if (reading_elements)
    {
      if (fill_bc_data)
      {
        long long int v[P4EST_CHILDREN / 2];
        long long int n;
        long long int el_num;
        int retval;

        if (num_elements >= *num_bc)
        {
          P4EST_LERROR("Encountered element that will not fit into"
                       " bc_to_vertex array. More elements than expected.\n");
          P4EST_FREE(line);
          return 1;
        }

        if (wrong_element)
        {
          retval = sscanf(line, "%lld, %lld, %lld, %lld", &el_num, &v[0], &v[1],
                          &v[2]);

          if (retval != 4)
          {
            P4EST_LERROR("Premature end of file");
            P4EST_FREE(line);
            return 1;
          }

          bc_el_label[num_elements] = el_num;
          for (n = 0; n < 3; ++n)
            bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + n] = v[n] - 1;
          bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + 3] = -1;
        }
        else
        {
          retval = sscanf(line, "%lld, %lld, %lld"
#ifdef P4_TO_P8
                                ", %lld, %lld"
#endif
                          ,
                          &el_num, &v[0], &v[1]
#ifdef P4_TO_P8
                          ,
                          &v[2], &v[3]
#endif
                          );

          if (retval != P4EST_CHILDREN / 2 + 1)
          {
            P4EST_LERROR("Premature end of file");
            P4EST_FREE(line);
            return 1;
          }

          bc_el_label[num_elements] = el_num;
          for (n = 0; n < P4EST_CHILDREN / 2; ++n)
            bc_to_vertex[P4EST_CHILDREN / 2 * num_elements + n] = v[n] - 1;
        }
      }
      ++num_elements;
    }
    else if (reading_elsets)
    {
      if (fill_bc_data)
      {

        char *p = line;
        while (*p)
        {
          if (isdigit(*p))
          {
            long val = strtol(p, &p, 10);
            element_label[ibc] = val;
            physical_type[ibc] = boundary_condition;
            ibc++;
          }
          else
          {
            p++;
          }
        }

        if (ibc > *num_bc)
        {
          P4EST_LERROR(
              "Encountered bc element that will not fit in bc array\n");
          P4EST_FREE(line);
          return 1;
        }
      }

      ++num_elsets;
    }
    else if (reading_elsets_volume)
    {
      if (fill_bc_data)
      {

        char *p = line;
        while (*p)
        {
          if (isdigit(*p))
          {
            long val = strtol(p, &p, 10);
            vc_el_label[ivc] = val;
            vc_el_type[ivc] = boundary_condition;
            ivc++;
          }
          else
          {
            p++;
          }
        }

        if (ivc > *num_elem)
        {
          P4EST_LERROR(
              "Encountered vc element that will not fit in vc array\n");
          P4EST_FREE(line);
          return 1;
        }
      }

      ++num_elsets_volume;
    }
    else if (reading_elements_volume)
    {
      if (fill_bc_data)
      {
        long long int v[P4EST_CHILDREN];
        long long int el_num;
        int retval;

        retval = sscanf(line, "%lld, %lld, %lld, %lld, %lld"
#ifdef P4_TO_P8
                              ", %lld, %lld, %lld, %lld"
#endif
                        ,
                        &el_num, &v[0], &v[1], &v[2], &v[3]
#ifdef P4_TO_P8
                        ,
                        &v[4], &v[5], &v[6], &v[7]
#endif
                        );

        if (retval != P4EST_CHILDREN + 1)
        {
          P4EST_LERROR("Premature end of file");
          P4EST_FREE(line);
          return 1;
        }

        vc_el_label_list[num_elements_volume] = el_num;
      }

      ++num_elements_volume;
    }

    ++lines_free;
    P4EST_FREE(line);
  }

  if (fill_bc_data)
  {
    for (int i = 0; i < *num_bc; i++)
    {
      for (int j = 0; j < *num_bc; j++)
      {
        if (element_label[i] == bc_el_label[j])
        {
          bc_physical_type[j] = physical_type[i];
          continue;
        }
      }
    }

    for (int i = 0; i < ivc; i++) // could be shorter loop
    {
      for (int j = 0; j < *num_elem; j++)
      {
        if (vc_el_label[i] == vc_el_label_list[j])
        {
          vc_el_label[i] = j;
          continue;
        }
      }
    }
  }
  *num_bc = num_elements;
  *num_elem = num_elements_volume;

  if (num_elsets == 0 || num_elements == 0)
  {
    P4EST_LERROR("No elements or bcs found in mesh file.\n");
    return -1;
  }
  else
  {
    return 0;
  }
}

// static int p4est_bc_read_inp (const char *filename, p4est_bc_info_t *bc_data)
static p4est_bc_info_t *p4est_bc_read_inp(const char *filename)
{
  FILE *fid = NULL;
  int retval;
  p4est_topidx_t num_bc = 0;
  p4est_topidx_t num_elem = 0;

  P4EST_GLOBAL_PRODUCTIONF("Reading BCs from %s\n", filename);

  fid = fopen(filename, "rb");
  if (fid == NULL)
  {
    P4EST_LERRORF("Failed to open %s\n", filename);
    goto dead;
  }

  rewind(fid);

  // read number of bc elements
  p4est_bc_read_inp_stream(fid, &num_bc, &num_elem, NULL, NULL, NULL, NULL,
                           NULL);

  rewind(fid);

  // allocate memory
  bc_data = malloc(sizeof(p4est_bc_info_t));
  bc_data->bc_physical_type = malloc(num_bc * sizeof(p4est_topidx_t));
  bc_data->bc_el_label = malloc(num_bc * sizeof(p4est_topidx_t));
  bc_data->bc_to_vertex =
      malloc(P4EST_CHILDREN / 2 * num_bc * sizeof(p4est_topidx_t));
  bc_data->num_bc = num_bc;
  bc_data->vc_el_type = malloc(num_elem * sizeof(p4est_topidx_t));
  bc_data->vc_el_label = malloc(num_elem * sizeof(p4est_topidx_t));

  memset(bc_data->vc_el_type, 0, sizeof(p4est_topidx_t) * num_elem);
  memset(bc_data->vc_el_label, 0, sizeof(p4est_topidx_t) * num_elem);

  // read boundary element data
  p4est_bc_read_inp_stream(fid, &num_bc, &num_elem, bc_data->bc_physical_type,
                           bc_data->bc_el_label, bc_data->bc_to_vertex,
                           bc_data->vc_el_type, bc_data->vc_el_label);

  int i = 0;
  for (i = 0; i < num_elem; i++)
  {
    if (bc_data->vc_el_type[i] == 0)
      break;
  }
  bc_data->nel_vc_sponge = i;

  // sort vertices
  for (int i = 0; i < num_bc; i++)
  {
    int nvert = P4EST_CHILDREN / 2;
    p4est_topidx_t v_bc[nvert];

    // get vertices from the boundary face
    for (int j = 0; j < nvert; j++)
    {
      v_bc[j] = bc_data->bc_to_vertex[nvert * i + j];
    }
    // sort the boundary vertices
    qsort(v_bc, nvert, sizeof(p4est_topidx_t), compare);

    // put them back in
    for (int j = 0; j < nvert; j++)
    {
      bc_data->bc_to_vertex[nvert * i + j] = v_bc[j];
    }
  }

  retval = fclose(fid);
  fid = NULL;
  if (retval)
  {
    P4EST_LERRORF("Failed to close %s\n", filename);
    goto dead;
  }

  return bc_data;

dead:
  /* clean up on error */
  if (fid != NULL)
  {
    fclose(fid);
  }
  return NULL;
}

// Function finds whether the element e is in the set s of length n.

static int is_in(p4est_topidx_t e, p4est_topidx_t *s, int n)
{
  long long int i;
  int match = 0;

  for (i = 0; i < n; i++)
  {
    if (e == s[i])
    {
      match = 1;
      break;
    }
  }
  return match;
}

// Function matches the vertices to the boundary element read from external
// meshfile
static int get_bc(p4est_topidx_t *v, p4est_bc_info_t *bc_data)
{
  long long int i;
  long long int nbf = bc_data->num_bc; // number of boundary faces
  int nvert = P4EST_CHILDREN / 2;
  p4est_topidx_t v_bc;
  int count = 0;
  int not_found_match = 1;

  // sort the input vertices
  qsort(v, nvert, sizeof(p4est_topidx_t), compare);

  // loop over all boundary faces extracted from the meshfile
  for (i = 0; i < nbf; i++)
  {

    count = 0;
    // get vertices from the boundary face (presorted)
    for (int j = 0; j < nvert; j++)
    {
      v_bc = bc_data->bc_to_vertex[nvert * i + j];
      if (v[j] == v_bc)
        count++;
    }
    if (count == nvert)
    {
      return bc_data->bc_physical_type[i];
      not_found_match = 0;
      break;
    }
  }
  if (not_found_match)
  {
    P4EST_LERROR("Did not find a match for external boundary condition! - "
                 "check .inp file for errors\n");
    return -1;
  }
}

// Function checks whether the element is in the sponge zone
static p4est_topidx_t get_vc(p4est_topidx_t qt, p4est_bc_info_t *bc_data)
{
  long long int i;
  int count = 0;
  int not_found_match = 1;

  for (i = 0; i < bc_data->nel_vc_sponge; i++)
  {
    if (qt == bc_data->vc_el_label[i])
    {
      return bc_data->vc_el_type[i];
      not_found_match = 0;
      break;
    }
  }
  if (not_found_match)
    return 0;
}

/* ------------------------------------------------------ */
/* THIS IS THE MAIN DRIVER FUNCTION CALLED FROM MOD_P4EST */
/* ------------------------------------------------------ */
void P2N_FUN(init)(int *is_cube, int *nnx, int *nny, int *nnz, int *nref_levs,
                   int *read_external_grid_flg, int *is_non_conforming_flg,
                   int *refine_elements, int *xperiodic_flg, int *yperiodic_flg,
                   int *zperiodic_flg, int *mFACE_LEN, int *mORIENT,
                   int *lrestoring_sponge)
{
  /* March 2015, Modified by Simone Marras to read external grids into p4est
     (abaqus format: see instructions on input file format in
     p4est_connectivity.h)
   */

  mpi_context_t *mpi = &mpi_context;
  p4est_refine_t refine_fn;
  p4est_t *p4est = NULL;
  p4est_connectivity_t *connectivity;
  p4est_quadrant_t *quad;
  int nx;
  int ny;
  int xperiodic, yperiodic;
  // int orient;
  int nref_levels;
  int read_external;
  p4est_locidx_t q;
  int is_geometry_cube;
  int is_non_conforming;
  int is_restoring;

  nx = *nnx;
  ny = *nny;
  xperiodic = *xperiodic_flg;
  yperiodic = *yperiodic_flg;
  FACE_LEN = *mFACE_LEN;
  ORIENT = *mORIENT;

  is_geometry_cube = *is_cube;

  nref_levels = *nref_levs;
  read_external = *read_external_grid_flg;
  is_non_conforming = *is_non_conforming_flg;
  is_restoring = *lrestoring_sponge;

  /* Select refinement function depending if non-conforming refinement allowed*/
  if (is_non_conforming)
  {
    refine_fn = refine_list_fn;
  }
  else
  {
    refine_fn = refine_uniform_fn;
  }

  /* Select p4est connectivity based on space_method*/
  p4est_connect_type_t connect_type;
  connect_type = P4EST_CONNECT_FULL;
  refine_level = nref_levels;

  /** Create initial grid */
  if (!stored_p4est)
  {

    /* Build or read grid */
    if (read_external == 0)
    {
#ifdef P4_TO_P8
      if (!is_geometry_cube)
      {
        // connectivity = p8est_connectivity_new_shell();
        connectivity = p8est_connectivity_new_cubed_sphere();
      }
      else
      {
        int zperiodic = *zperiodic_flg;
        int nz;
        nz = *nnz;
        connectivity = p4est_connectivity_new_brick(nx, ny, nz, xperiodic,
                                                    yperiodic, zperiodic);
      }
#else
      p4est_connectivity_t *conn_init;
      if (!is_geometry_cube)
      {
        conn_init = p4est_connectivity_new_cubed_sphere();
        connectivity = p4est_connectivity_refine(conn_init, nx);
        p4est_connectivity_destroy(conn_init);
      }
      else
      {
        connectivity =
            p4est_connectivity_new_brick(nx, ny, xperiodic, yperiodic);
      }
#endif
    }
    else
    {
      connectivity = p4est_connectivity_read_inp("EXTERNAL_MESH.inp");
      p4est_connectivity_memory_used(connectivity);

#ifdef P4EST_WITH_METIS
//      p4est_connectivity_reorder (mpi->mpicomm, 0, connectivity,
//      P4EST_CONNECT_FACE);
#endif
      /* P4EST_WITH_METIS */

      bc_data = p4est_bc_read_inp("EXTERNAL_MESH.inp");
    }

    /* Create a forest that is not refined; it consists of the root octant.
     * The p4est_new_ext function can take a startlevel for a load-balanced
     * initial uniform refinement.  Here we refine adaptively instead. */
    p4est = (p4est_t *)p4est_new_ext(mpi->mpicomm, connectivity, 0, nref_levels,
                                     1, sizeof(user_data_t), init_fn, NULL);
  }
  else
  {
    p4est = stored_p4est;
    connectivity = p4est->connectivity;
  }

  /* Set refinement refinement flag */
  if (is_non_conforming)
  {
    /*Fill refine flag in user_data*/
    {
      user_data_t *data;

      /* Use do_refine array to set quadrant user_data, which holds refinemenr
       * flags  */
      for (q = 0; q < p4est->local_num_quadrants; ++q)
      {
        /*  find which quadrant */
        quad = p4est_mesh_quadrant_cumulative(p4est, q, 0, NULL);
        data = (user_data_t *)quad->p.user_data;
        /*  assign user_int from do_refine */
        data->iref = refine_elements[q];
      }
    }
  }

  /* refine, balance and partition*/
  {
    p4est_refine_ext(p4est, 1, -1, refine_fn, NULL, quadrant_init_replace);
    p4est_balance(p4est, connect_type, init_fn);
    p4est_partition(p4est, 1, NULL);
  }

/* Write vtk file */
#ifdef VTK_OUTPUT
  p4est_vtk_write_file(p4est, NULL, "p4est_test_00");
#endif

  stored_p4est = p4est;
}

void P2N_FUN(fill_data)(int *nnop, double *xgl, int *nis_dg, int *iboundary,
                        int *read_external_grid_flg, p4esttonuma_t **p2n,
                        int *lrestoring_sponge)
{
  mpi_context_t *mpi = &mpi_context;
  p4est_t *p4est = stored_p4est;
  p4est_ghost_t *ghost;
  p4est_lnodes_t *lnodes, *fnodes, *cnodes;
  p4est_mesh_t *mesh;
  p4est_locidx_t num_physical_boundary_faces;
  p4est_locidx_t num_processor_boundary_faces;
  p4est_quadrant_t *quad;
  int nop = *nnop;
  size_t zz;
  int i;
  int num_send_recv_total, num_send_recv_offset, num_nbh;
  p4est_locidx_t ll, q, nq, sk, skb;
  int f, nf, j, nv, h;
  int is_dg = *nis_dg;
  int read_external = *read_external_grid_flg;

  int is_restoring = *lrestoring_sponge;

  /* Create the ghost layer to learn about parallel neighbors. */
  p4est_connect_type_t connect_type;
  connect_type = P4EST_CONNECT_FULL;
  ghost = p4est_ghost_new(p4est, connect_type);

  /* Create face numbering for DG and CDG method */
  fnodes = p4est_lnodes_new(p4est, ghost, -1);

  /* Create a node numbering for continuous linear finite elements. */
  lnodes = p4est_lnodes_new(p4est, ghost, nop);

  /* Create a new balanced mesh with quad-to-tree structure*/
  mesh = p4est_mesh_new_ext(p4est, ghost, 1, 0, connect_type);

  /* allocate p2n */
  *p2n = malloc(sizeof(p4esttonuma_t));

  /* nop: order of polynomial in 1D */
  (*p2n)->nop = nop;

  /* npoin: number of grid points on the processor */
  (*p2n)->npoin = lnodes->num_local_nodes;

  /* nfaces: number of faces on the processor */
  (*p2n)->nface = fnodes->num_local_nodes;

  /* nelem: number of elements on the processor */
  (*p2n)->nelem = p4est->local_num_quadrants;

/* npoin_dg: number of dg points */
#ifdef P4_TO_P8
  (*p2n)->npoin_dg = (*p2n)->nelem * (nop + 1) * (nop + 1) * (nop + 1);
#else
  (*p2n)->npoin_dg = (*p2n)->nelem * (nop + 1) * (nop + 1);
#endif

  /* nface_dg: number of dg faces */
  (*p2n)->nface_dg = P4EST_FACES * p4est->local_num_quadrants;

  if (is_restoring)
  {
    (*p2n)->vc_el_type = malloc(sizeof(p4est_topidx_t) * (*p2n)->nelem);
    memset((*p2n)->vc_el_type, 0, sizeof(p4est_topidx_t) * (*p2n)->nelem);
  }

  /* Get coordinates */
  fill_coordinates(nop + 1, xgl, p4est, lnodes, *p2n);

  /* Get intma_table *
   * intma_table(1:nglx,1:ngly,1:nglz,1:nelem): gives point number (1...npoin)
   * for  *
   * each point in each element (ngl*=nop*+1, nop: polynomial degree)         */
  (*p2n)->intma_table = malloc(sizeof(p4est_locidx_t) * (*p2n)->npoin_dg);
  for (ll = 0; ll < (*p2n)->npoin_dg; ++ll)
    (*p2n)->intma_table[ll] = lnodes->element_nodes[ll] + 1;

  /* Count the number of elements that touch non-conforming elements */
  p4est_locidx_t *EToNC = (*p2n)->EToNC =
      malloc(sizeof(p4est_locidx_t) * p4est->local_num_quadrants);
  p4est_locidx_t nNC = 0;
  for (int n = 0; n < p4est->local_num_quadrants; ++n)
    if (lnodes->face_code[n])
    {
      EToNC[n] = nNC;
      ++nNC;
    }
    else
      EToNC[n] = -1;
  (*p2n)->nNC = nNC;

  /* For NC elements, determine how the faces/edges are hanging */
  int *NC_face = (*p2n)->NC_face = malloc(sizeof(int) * nNC * P4EST_FACES);
#ifdef P4_TO_P8
  int *NC_edge = (*p2n)->NC_edge = malloc(sizeof(int) * nNC * P8EST_EDGES);
#else
  (*p2n)->NC_edge = NULL;
#endif

  int m = 0;
  for (int n = 0; n < p4est->local_num_quadrants; ++n)
    if (lnodes->face_code[n])
    {
      p4est_lnodes_decode(lnodes->face_code[n], NC_face + m * P4EST_FACES
#ifdef P4_TO_P8
                          ,
                          NC_edge + m * P8EST_EDGES
#endif
                          );
      ++m;
    }

  if (m != nNC)
    SC_ABORT("problem with counting face_code");

  /**
   * Face id
   */
  p4est_locidx_t *intma_face =
      malloc(sizeof(p4est_locidx_t) * (*p2n)->nface_dg);
  for (ll = 0; ll < (*p2n)->nface_dg; ++ll)
    intma_face[ll] = fnodes->element_nodes[ll] + 1;

  /* Switch between nodes/faces enumeration *
   * for CG and DG respectively             */
  cnodes = is_dg ? fnodes : lnodes;

  num_nbh = (int)cnodes->sharers->elem_count;

  /*malloc processor arrays*/
  (*p2n)->num_send_recv = malloc(sizeof(int) * num_nbh);
  (*p2n)->nbh_proc = malloc(sizeof(int) * num_nbh);

  /*
   * Fill the sharing information for CG and DG
   */
  for (zz = 0, num_send_recv_total = 0, num_nbh = 0;
       zz < cnodes->sharers->elem_count; ++zz)
  {
    p4est_lnodes_rank_t *lrank =
        p4est_lnodes_rank_array_index(cnodes->sharers, zz);

    if (lrank->rank != mpi->mpirank)
    {
      (*p2n)->nbh_proc[num_nbh] = lrank->rank + 1;
      (*p2n)->num_send_recv[num_nbh] = (int)lrank->shared_nodes.elem_count;

      num_send_recv_total += (int)lrank->shared_nodes.elem_count;
      num_nbh += 1;
    }
  }

  /*
   * num_send_recv_total: total number of points (or faces) to be communicated
   * with all
   * neighboring processors
   */
  (*p2n)->num_send_recv_total = num_send_recv_total;

  /* num_nbh: number of neighboring processors */
  (*p2n)->num_nbh = num_nbh;

  /* nbh_send_recv(1:num_send_recv_total): lcoal number of the point that must
   * be sent/received sorted in a global ordering (e.g.:
   * nbh_send_recv(1:num_send_recv(1)) are the points to be communicated with
   * processor nbh_proc(1).
   * nbh_send_recv(num_send_recv(1)+1:num_send_recv(1)+num_send_recv(2)) are
   * the points to be communicated with processor nbh_proc(2)...) */

  (*p2n)->nbh_send_recv = malloc(sizeof(p4est_locidx_t) * num_send_recv_total);
  (*p2n)->nbh_send_recv_multi =
      malloc(sizeof(p4est_locidx_t) * num_send_recv_total);
  (*p2n)->nbh_send_recv_half =
      malloc(sizeof(p4est_locidx_t) * num_send_recv_total * P4EST_HALF);

  memset((*p2n)->nbh_send_recv_multi, 0,
         sizeof(p4est_locidx_t) * num_send_recv_total);
  memset((*p2n)->nbh_send_recv_half, 0,
         sizeof(p4est_locidx_t) * num_send_recv_total * P4EST_HALF);

  (*p2n)->face_type = malloc(sizeof(p4est_locidx_t) * (*p2n)->nface);

  memset((*p2n)->face_type, 0, sizeof(p4est_locidx_t) * (*p2n)->nface);

  for (zz = 0, num_send_recv_offset = 0; zz < cnodes->sharers->elem_count; ++zz)
  {
    p4est_lnodes_rank_t *lrank =
        p4est_lnodes_rank_array_index(cnodes->sharers, zz);
    const int nshared = (int)lrank->shared_nodes.elem_count;
    int jn;

    if (lrank->rank != mpi->mpirank)
    {
      for (jn = 0; jn < nshared; ++jn)
      {
        int gl = (int)*(p4est_locidx_t *)sc_array_index_int(
            &lrank->shared_nodes, jn);
        (*p2n)->nbh_send_recv[num_send_recv_offset + jn] = gl + 1;
      }
      num_send_recv_offset += nshared;
    }
  }

  /*/
    printf("NumSend %d\n",(*p2n)->num_send_recv_total);
    for(i = 0; i < (*p2n)->num_send_recv_total; ++i)
    printf("%d ",(*p2n)->nbh_send_recv[i]);
    printf("\n");
  */

  /* Loop over all processor boundary faces and mark them accordingly
   * (face_type=2)*/
  for (zz = 0; zz < fnodes->sharers->elem_count; ++zz)
  {
    p4est_lnodes_rank_t *frank =
        p4est_lnodes_rank_array_index(fnodes->sharers, zz);
    const int nshared = (int)frank->shared_nodes.elem_count;
    int jn;

    if (frank->rank != mpi->mpirank)
    {
      for (jn = 0; jn < nshared; ++jn)
      {
        int gl = (int)*(p4est_locidx_t *)sc_array_index_int(
            &frank->shared_nodes, jn);
        (*p2n)->face_type[gl] = 2;
      }
    }
  }

  /*/
    for(zz = 0 ; zz < num_nbh; zz++) {
    printf("[p4estnuma] Processor %d: neighbor %d total send recv %d\n",
  mpi->mpirank,
    (*p2n)->nbh_proc[zz], (*p2n)->num_send_recv[zz]);
    }

    for(zz = 0 ; zz < num_send_recv_total; zz++) {
    printf("[p4estnuma] Processor %d: zz %d, nbh_send_recv[zz]
  %d\n",mpi->mpirank,zz,(*p2n)->nbh_send_recv[zz]);
    }
  */

  /* printf("NumVertices %d\n",p4est->connectivity->num_vertices);
  for(int i = 0; i < p4est->connectivity->num_vertices; ++i)
    {
      printf("%d\t",i);
      for(int j = 0; j<3; ++j) printf("%f
  \t",p4est->connectivity->vertices[3*i+j]);
      printf("\n");
    }
  printf("tree_to_vertex\n");
  for(int i = 0; i < p4est->connectivity->num_trees; ++i)
    {
      printf("%d\t",i);
      for(int j = 0; j<4; ++j) printf("%d
  \t",p4est->connectivity->tree_to_vertex[4*i+j]);
      printf("\n");
    }
  printf("NumTrees %d\n",p4est->connectivity->num_trees);
  printf("tree_to_face\n");
  for(int i = 0; i < p4est->connectivity->num_trees; ++i)
    {
      printf("%d\t",i);
      for(int j = 0; j<4; ++j) printf("%d
  \t",p4est->connectivity->tree_to_face[4*i+j]);
      printf("\n");
    }
  printf("tree_to_tree\n");
  for(int i = 0; i < p4est->connectivity->num_trees; ++i)
    {
      printf("%d\t",i);
      for(int j = 0; j<4; ++j) printf("%d
  \t",p4est->connectivity->tree_to_tree[4*i+j]);
      printf("\n");
    }
  */

  /*
   * Count boundary faces
   */
  num_physical_boundary_faces = 0;
  // TODO: Understand why this works for CG but not the routine below
  num_processor_boundary_faces = 0;
  for (q = 0; q < mesh->local_num_quadrants; ++q)
  {
    for (f = 0; f < P4EST_FACES; ++f)
    {
      nq = mesh->quad_to_quad[P4EST_FACES * q + f];
      nv = mesh->quad_to_face[P4EST_FACES * q + f];
      nf = nv % P4EST_FACES;

      // If I am my neighbor then boundary
      if (nq == q && nf == f)
      {
        ++num_physical_boundary_faces;
      }
      // If my neighbors face is negative I'm hanging
      else if (nv < 0)
      {
        p4est_locidx_t *nqs;
        int h;
        nqs = sc_array_index(mesh->quad_to_half, nq);
        // Check whether neighbor children are local or not
        for (h = 0; h < P4EST_HALF; ++h)
        {
          if (nqs[h] >= mesh->local_num_quadrants)
          {
            ++num_processor_boundary_faces;
          }
        }
      }
      // I'm conforming or a little guy, so see if my neighbor is local or not
      // (this count is different than the one done below, and thus a bug)
      else if (nq >= mesh->local_num_quadrants)
      {
        ++num_processor_boundary_faces;
      }
    }
  }

  /*
   * nbsido: number of (physical) domain boundary faces on the processor
   */
  (*p2n)->nbsido = num_physical_boundary_faces;

  /*
   * nboun: number of domain and processor boundary faces on the processor
   * (processor boundary faces + domain boundary faces, see
   * domain_decomp_metis.f90, line 1440) */

  (*p2n)->nboun =
      /*num_physical_boundary_faces + */ num_processor_boundary_faces;

/*face orientation transformation array p4est->numa*/
#ifdef P4_TO_P8
  static const int transform[P4EST_FACES] = {4, 5, 2, 3, 0, 1};
#else
  static int transform[P4EST_FACES];
  // yz
  if (ORIENT == 0)
  {
    transform[0] = 2;
    transform[1] = 3;
    transform[2] = 0;
    transform[3] = 1;
    // xz
  }
  else if (ORIENT == 1)
  {
    transform[0] = 4;
    transform[1] = 5;
    transform[2] = 0;
    transform[3] = 1;
    // xy
  }
  else
  {
    transform[0] = 4;
    transform[1] = 5;
    transform[2] = 2;
    transform[3] = 3;
  }
#endif

  /*
   * bsido(1:6,1:nbsido): domain boundary data:
   * bsido(1:4,i): point number (1...npoin) of the four corners of the face i
   *               (counter-clock wise when seen from outside of the domain)
   * bsido(5,i): element (1...nelem) to which face i belongs (only one because
   *             boundary face!)
   * bsido(6,i): boundary condition flag for face i (=4 for solid wall boundary)
   */

  (*p2n)->bsido =
      malloc(sizeof(p4est_locidx_t) * 6 * num_physical_boundary_faces);
  memset((*p2n)->bsido, 0, sizeof(p4est_locidx_t) * 6 * (*p2n)->nbsido);

  /*
   * face(1:8, 1:nface): face data:
   * face(1:4,i): point number (1...npoin) of the four corners of the face i
   *               (counter-clock wise when seen from outside of the domain)
   * face(5,i): local face relative to left element
   * face(6,i): local face relative to right element
   * face(7,i): element to left
   * face(8,i): element to the right or boundary condition flag (-bsido(6,sk))
   *
   * for non-conforming faces entries 8-11 store four half-size neighbor
   * elements */

  (*p2n)->face = malloc(sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);
  memset((*p2n)->face, 0, sizeof(p4est_locidx_t) * FACE_LEN * (*p2n)->nface);

#ifdef P4_TO_P8
  static const int CONFORMING_LIM = 24;
#else
  static const int CONFORMING_LIM = 8;
#endif

  int count_matching_faces = 0;
  int count_mismatched_bcs = 0;

  num_processor_boundary_faces = 0;
  for (q = 0, skb = 0; q < mesh->local_num_quadrants; ++q)
  {
    for (f = 0; f < P4EST_FACES; ++f)
    {
      nq = mesh->quad_to_quad[P4EST_FACES * q + f];
      nv = mesh->quad_to_face[P4EST_FACES * q + f];
      nf = nv % P4EST_FACES;
      sk = intma_face[P4EST_FACES * q + f] - 1;

      if (nq == q && nf == f)
      { // physical boundary faces

        // (*p2n)->face_type[sk] = 4;

        p4est_topidx_t v[P4EST_CHILDREN / 2];

        // Get the vertices of a matching tree face

        if (read_external)
        {
          for (i = 0; i < P4EST_CHILDREN / 2; i++)
          {
            int qtt = mesh->quad_to_tree[q];
            int vv = p4est_face_corners[f][i];
            v[i] =
                p4est->connectivity->tree_to_vertex[P4EST_CHILDREN * qtt + vv];
          }

          // check whether the vertices match any of the bc elements

          (*p2n)->face_type[sk] = get_bc(v, bc_data);
          if ((*p2n)->face_type[sk] > 0)
          {
            count_matching_faces++;
          }
          else if ((*p2n)->face_type[sk] == -1)
          {
#ifdef P4_TO_P8
            printf("Missed BC: %d %d %d %d\n", v[0], v[1], v[2], v[3]);
#else
            printf("Missed BC: %d %d\n", v[0], v[1]);
#endif
            //	      printf("Coords: (%f %f %f),(%f %f %f),(%f %f %f),(%f %f
            //%f)\n",
            //     (*p2n)->coord_dg
            count_mismatched_bcs++;
          }

          (*p2n)->face[FACE_LEN * sk + 7] = -(*p2n)->face_type[sk];
          (*p2n)->bsido[6 * skb + 5] = -(*p2n)->face_type[sk];
        }
        else
        {
          (*p2n)->face[FACE_LEN * sk + 7] = -iboundary[transform[f]];
          (*p2n)->bsido[6 * skb + 5] = iboundary[transform[f]];
        }

        (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
        (*p2n)->face[FACE_LEN * sk + 5] = 0;
        (*p2n)->face[FACE_LEN * sk + 6] = q + 1;

        (*p2n)->bsido[6 * skb + 4] = q + 1;
        //

        ++skb;
      }
      else if ((nv >= 0) && (nv < CONFORMING_LIM))
      { // conforming face
        if (nq >= mesh->local_num_quadrants)
        { // processor boundary face

          (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
          (*p2n)->face[FACE_LEN * sk + 5] = 0;
          (*p2n)->face[FACE_LEN * sk + 6] = q + 1;
          (*p2n)->face[FACE_LEN * sk + 7] =
              0; // nq - mesh->local_num_quadrants +1; //HACK
          (*p2n)->face_type[sk] = 2;

          int nbh = mesh->ghost_to_proc[nq - mesh->local_num_quadrants];
          int offset = 0;

          for (i = 0; i < num_nbh; i++)
          { // go over all neighbor procs
            if ((*p2n)->nbh_proc[i] == nbh + 1)
            { // find neighbor on the list
              // go over faces to be send to this proc
              for (j = offset; j < offset + (*p2n)->num_send_recv[i]; j++)
              {
                if ((*p2n)->nbh_send_recv[j] == sk + 1)
                {
                  // increment the multiplicity of this face
                  (*p2n)->nbh_send_recv_multi[j]++;
                  ++num_processor_boundary_faces;
                }
              }
            }
            offset += (*p2n)->num_send_recv[i];
          }
        }
        else if (q < nq || (q == nq && f < nf))
        { // left element
          (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
          (*p2n)->face[FACE_LEN * sk + 6] = q + 1;
          (*p2n)->face_type[sk] = 1;
        }
        else
        { // right element
          (*p2n)->face[FACE_LEN * sk + 5] = transform[f] + 1;
          (*p2n)->face[FACE_LEN * sk + 7] = q + 1;
          (*p2n)->face_type[sk] = 1;
        }
      }
      else if (nv < 0)
      { // non-conforming face: half-size neighbors
        (*p2n)->face[FACE_LEN * sk + 4] = transform[f] + 1;
        (*p2n)->face[FACE_LEN * sk + 6] = q + 1; // i am the parent

        if ((*p2n)->face_type[sk] == 2)
          (*p2n)->face_type[sk] = 21; // mark as non-conforming processor edge
                                      //         else
        //           (*p2n)->face_type[sk] = 1; //internal

        p4est_locidx_t *nqs;
        int h;
        nqs = sc_array_index(mesh->quad_to_half, nq);

        for (h = 0; h < P4EST_HALF; ++h)
          if (nqs[h] >= mesh->local_num_quadrants)
          {
            //    (*p2n)->face_type[sk] = 21; //mark as non-conforming processor
            //    edge
            int nbh = mesh->ghost_to_proc[nqs[h] - mesh->local_num_quadrants];
            int offset = 0;

            for (i = 0; i < num_nbh; i++)
            { // go over all neighbor procs

              if ((*p2n)->nbh_proc[i] == nbh + 1)
              { // find neighbor on the list
                // go over faces to be send to this proc
                for (j = offset; j < offset + (*p2n)->num_send_recv[i]; j++)
                {
                  if ((*p2n)->nbh_send_recv[j] == sk + 1)
                  {
                    // increment the multiplicity of this face
                    (*p2n)->nbh_send_recv_multi[j]++;
                    ++num_processor_boundary_faces;
                    (*p2n)->nbh_send_recv_half[j * P4EST_HALF + h] = 1;
                  }
                }
              }

              offset += (*p2n)->num_send_recv[i];
            }
          }
      }
      else if (nv >= CONFORMING_LIM)
      { // non-conforming face: real-size neighbor
        (*p2n)->face[FACE_LEN * sk + 5] = transform[f] + 1;
        h = nv / CONFORMING_LIM - 1; // find my sub-face number
        (*p2n)->face[FACE_LEN * sk + 7 + h] =
            q + 1; // put me in the right child space
        if ((*p2n)->face_type[sk] == 2)
          (*p2n)->face_type[sk] = 21;
        //         else
        //           (*p2n)->face_type[sk] = 1; //internal

        // If this quadrant is a ghost then set up some comm info
        if (nq >= mesh->local_num_quadrants)
        {
          int nbh = mesh->ghost_to_proc[nq - mesh->local_num_quadrants];
          int offset = 0;

          for (i = 0; i < num_nbh; i++)
          { // go over all neighbor procs
            if ((*p2n)->nbh_proc[i] == nbh + 1)
            { // find neighbor on the list
              //(*p2n)->face_type[sk] = 21;
              // go over faces to be send to this proc
              for (j = offset; j < offset + (*p2n)->num_send_recv[i]; j++)
              {
                if ((*p2n)->nbh_send_recv[j] == sk + 1)
                {
                  // increment the multiplicity of this face
                  (*p2n)->nbh_send_recv_multi[j]++;
                  ++num_processor_boundary_faces;
                }
              }
            }
            offset += (*p2n)->num_send_recv[i];
          }
        }
      }
    }

    // Check whether this quadrant is in the sponge zone
    // to do that I need to check whether tree vertices of local elements
    // match the vertices of sponge elements in the external mesh
    if (read_external && is_restoring)
    {
      int qtt = mesh->quad_to_tree[q];
      (*p2n)->vc_el_type[q] = get_vc(qtt, bc_data);
    }
  }

  if (count_mismatched_bcs > 0)
    printf("Found %d mismatched BCs and %d matched.\n", count_mismatched_bcs,
           count_matching_faces);

  if (is_dg)
    (*p2n)->nboun = num_processor_boundary_faces;

  // mark parent side of nc ip face
  for (sk = 0; sk < (*p2n)->nface; ++sk)
  {
    if (((*p2n)->face_type[sk] == 21) && ((*p2n)->face[FACE_LEN * sk + 6] > 0))
    {
      (*p2n)->face_type[sk] = 12;
    }
  }

  /* add the missing subfaces from nc faces from parent side */
  //     (*p2n)->nboun += num_processor_boundary_subfaces;

  /*   //faces */
  /*   for(sk = 0; sk < (*p2n)->nface; ++sk) */
  /*   { */
  /*       printf("Proc. %4d Face %4d. ", mpi->mpirank, sk+1); */
  /*       for(f = 0;f < FACE_LEN;f++) */
  /*         printf("%7d ",(*p2n)->face[sk*FACE_LEN+f]); */
  /*       printf("%7d ", (*p2n)->face_type[sk]); */
  /*       printf("\n"); */
  /*   } */
  /*   //bsido */
  /*   for(sk = 0; sk < (*p2n)->nbsido; ++sk) */
  /*   { */
  /*       printf("%4d. ", sk+1); */
  /*       for(f = 0;f < 6;f++) */
  /*         printf("%7d ",(*p2n)->bsido[sk*6+f]); */
  /*       printf("\n"); */
  /*   } */

  // Save a list of sponge elements

  if (read_external)
  {
    free(bc_data->vc_el_label);
    free(bc_data->vc_el_type);
    free(bc_data->bc_el_label);
    free(bc_data->bc_physical_type);
    free(bc_data->bc_to_vertex);
    free(bc_data);
  }

  // Following is for multirate only
  (*p2n)->level_list =
      malloc(sizeof(p4est_locidx_t) * mesh->local_num_quadrants);
  (*p2n)->plist = malloc(sizeof(p4est_locidx_t) * mesh->local_num_quadrants);

  int *temp_list = malloc(sizeof(p4est_locidx_t) * mesh->local_num_quadrants);

  memset((*p2n)->plist, 0, sizeof(p4est_locidx_t) * mesh->local_num_quadrants);

  int neighbor = 0;
  for (q = 0; q < mesh->local_num_quadrants; ++q)
  {
    quad = p4est_mesh_quadrant_cumulative(p4est, q, 0, NULL);
    // save quad level
    (*p2n)->level_list[q] = quad->level;
    // if refined add to partition list
    if (quad->level > 0)
    {
      (*p2n)->plist[q] = 1;
      for (f = 0; f < P4EST_FACES; ++f)
      {
        sk = intma_face[P4EST_FACES * q + f] - 1;
        for (i = 0; i < FACE_LEN - 6; i++)
        {
          neighbor = (*p2n)->face[FACE_LEN * sk + 6 + i];
          if (neighbor > 0)
            (*p2n)->plist[neighbor - 1] = 1;
        }
      }
    }
  }

  /* add neighbors to partition list again to get two layer buffer and corner
   * elements */
  for (q = 0; q < mesh->local_num_quadrants; ++q)
  {
    temp_list[q] = (*p2n)->plist[q];
  }

  for (q = 0; q < mesh->local_num_quadrants; ++q)
  {
    if (temp_list[q] > 0)
    {
      for (f = 0; f < P4EST_FACES; ++f)
      {
        sk = intma_face[P4EST_FACES * q + f] - 1;
        for (i = 0; i < FACE_LEN - 6; i++)
        {
          neighbor = (*p2n)->face[FACE_LEN * sk + 6 + i];
          if (neighbor > 0)
            (*p2n)->plist[neighbor - 1] = 1;
        }
      }
    }
  }

  /*free intma_face*/
  free(temp_list);
  free(intma_face);

  /* Destroy the ghost structure -- no longer needed after node creation. */
  p4est_mesh_destroy(mesh);
  mesh = NULL;

  p4est_ghost_destroy(ghost);
  ghost = NULL;

  p4est_lnodes_destroy(fnodes);
  fnodes = NULL;

  p4est_lnodes_destroy(lnodes);
  lnodes = NULL;

  stored_p4est = p4est;
}

/* Function passes to NUMA some scalars describing the mesh, which are used for
   allocating arrays.
   It is called from MOD_P4EST */
void P2N_FUN(get_mesh_scalars)(p4esttonuma_t **p2n, int *npoin, int *nelem,
                               int *num_nbh, int *num_send_recv_total,
                               int *nbsido, int *nface, int *nboun, int *nNC)
{
  *npoin = (*p2n)->npoin;     /* Total number of CG nodes (local * non-local) */
  *nelem = (*p2n)->nelem;     /* Number of local elements */
  *num_nbh = (*p2n)->num_nbh; /* Number of neighboring ranks */
  *nbsido = (*p2n)->nbsido;   /* Number of faces on the physical boundary */
  *nface = (*p2n)->nface;     /* Number of unique CG faces on the processor */
  *nboun = (*p2n)->nboun;     /* Number of small faces (NC counted separately)*/
  /* Number of CG faces shared with neighboring ranks (NC counted together) */
  *num_send_recv_total = (*p2n)->num_send_recv_total;
  *nNC = (*p2n)->nNC;
}

/* Function passes to NUMA array data describing mesh. Called from MOD_P4EST*/
void P2N_FUN(get_mesh_arrays)(p4esttonuma_t **p2n, real *coord,
                              int *intma_table, int *NC_face, int *NC_edge,
                              int *EToNC, int *face, int *face_type,
                              int *num_send_recv, int *nbh_proc,
                              int *nbh_send_recv, int *nbh_send_recv_multi,
                              int *nbh_send_recv_half, int *bsido,
                              int *level_list, int *plist, int *nis_cgc,
                              int *vc_el_type, int *lrestoring_sponge)
{
  int i;
  int is_cgc;
  int is_restoring = *lrestoring_sponge;
  is_cgc = *nis_cgc;

  if (!is_cgc)
  {
    for (i = 0; i < 3 * (*p2n)->npoin_dg; ++i)
      coord[i] = (*p2n)->coord_dg[i];
  }
  else
  {
    for (i = 0; i < 3 * (*p2n)->npoin; ++i)
      coord[i] = (*p2n)->coord_cg[i];
  }

  for (i = 0; i < (*p2n)->npoin_dg; ++i)
    intma_table[i] = (*p2n)->intma_table[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    num_send_recv[i] = (*p2n)->num_send_recv[i];

  for (i = 0; i < (*p2n)->num_nbh; ++i)
    nbh_proc[i] = (*p2n)->nbh_proc[i];

  for (i = 0; i < (*p2n)->num_send_recv_total; ++i)
    nbh_send_recv[i] = (*p2n)->nbh_send_recv[i];

  for (i = 0; i < (*p2n)->num_send_recv_total; ++i)
    nbh_send_recv_multi[i] = (*p2n)->nbh_send_recv_multi[i];

  for (i = 0; i < (*p2n)->num_send_recv_total * P4EST_HALF; ++i)
    nbh_send_recv_half[i] = (*p2n)->nbh_send_recv_half[i];

  for (i = 0; i < 6 * (*p2n)->nbsido; ++i)
    bsido[i] = (*p2n)->bsido[i];

  for (i = 0; i < FACE_LEN * (*p2n)->nface; ++i)
    face[i] = (*p2n)->face[i];

  for (i = 0; i < (*p2n)->nface; ++i)
    face_type[i] = (*p2n)->face_type[i];

  for (i = 0; i < (*p2n)->nelem; ++i)
    level_list[i] = (*p2n)->level_list[i];

  for (i = 0; i < (*p2n)->nelem; ++i)
    plist[i] = (*p2n)->plist[i];

  for (i = 0; i < (*p2n)->nelem; ++i)
    EToNC[i] = (*p2n)->EToNC[i] + 1;

  for (i = 0; i < (*p2n)->nNC * P4EST_FACES; ++i)
    NC_face[i] = (*p2n)->NC_face[i];

#ifdef P4_TO_P8
  for (i = 0; i < (*p2n)->nNC * P8EST_EDGES; ++i)
    NC_edge[i] = (*p2n)->NC_edge[i];
#endif

  if (is_restoring)
  {
    for (i = 0; i < (*p2n)->nelem; ++i)
      vc_el_type[i] = (*p2n)->vc_el_type[i];
  }
}

/* Function frees memory used by P4EST*/
void P2N_FUN(free)(p4esttonuma_t **p2n, int *lrestoring_sponge)
{

  int is_restoring = *lrestoring_sponge;

  if (is_restoring)
    free((*p2n)->vc_el_type);
  free((*p2n)->EToNC);
  free((*p2n)->NC_face);
  if ((*p2n)->NC_edge)
    free((*p2n)->NC_edge);
  free((*p2n)->coord_dg);
  free((*p2n)->coord_cg);
  free((*p2n)->intma_table);
  free((*p2n)->num_send_recv);
  free((*p2n)->nbh_proc);
  free((*p2n)->nbh_send_recv);
  free((*p2n)->nbh_send_recv_multi);
  free((*p2n)->nbh_send_recv_half);
  free((*p2n)->bsido);
  free((*p2n)->face);
  free((*p2n)->face_type);
  free((*p2n)->level_list);
  free((*p2n)->plist);
  free(*p2n);
}

/* START OF BFAM CODE */
/*
 * The following code is extracted from/based on the bfam project:
 *
 *    https://github.com/bfam/bfam
 *
 * Authors:
 *   Lucas C Wilcox
 *   Jeremy E Kozdon
 */

static void mark_neighbors_volume_fun(p4est_iter_volume_info_t *info,
                                      void *user_data)
{
  user_data_t *qdata = info->quad->p.user_data;
  qdata->oref = qdata->iref;
}

static void mark_neighbors_corner_fun(p4est_iter_corner_info_t *info,
                                      void *user_data)
{
  user_data_t *ghost_data = (user_data_t *)user_data;
  sc_array_t *sides = &(info->sides);
  p4est_iter_corner_side_t *side;
  int iref = 0;
  int oref = 0;
  int8_t lvl = 0;

  for (int k = 0; k < sides->elem_count; ++k)
  {
    side = p4est_iter_cside_array_index_int(sides, k);
    if (side->is_ghost)
      oref = ghost_data[side->quadid].iref;
    else
      oref = ((user_data_t *)side->quad->p.user_data)->oref;

    if (oref > 0)
    {
      iref = 1;
      lvl = side->quad->level > lvl ? side->quad->level : lvl;
    }
  }

  if (iref > 0)
    for (int k = 0; k < sides->elem_count; ++k)
    {
      side = p4est_iter_cside_array_index_int(sides, k);
      if (side->is_ghost)
        continue;
      if (lvl >= side->quad->level)
      {
        user_data_t *qdata = side->quad->p.user_data;
        qdata->iref = 1;
      }
    }
}

void P2N_FUN(mark_neighbors_p4est)(int32_t *num_iterations)
{
  p4est_t *p4est = stored_p4est;
  p4est_ghost_t *ghost;
  p4est_connect_type_t connect_type;
  connect_type = P4EST_CONNECT_FULL;
  ghost = p4est_ghost_new(p4est, connect_type);
  user_data_t *ghost_data = NULL;

  ghost_data =
      (user_data_t *)malloc(sizeof(user_data_t) * ghost->ghosts.elem_count);

  for (int k = 0; k < *num_iterations; ++k)
  {
    p4est_ghost_exchange_data(p4est, ghost, ghost_data);
    p4est_iterate(p4est, ghost, (void *)ghost_data, mark_neighbors_volume_fun,
                  NULL,
#ifdef P4_TO_P8
                  NULL,
#endif
                  mark_neighbors_corner_fun);
  }

  free(ghost_data);

  p4est_ghost_destroy(ghost);
  ghost = NULL;
}

static void mark_element_fun(p4est_iter_volume_info_t *info, void *user_data)
{
  void **data = user_data;

  int *k = (int *)data[0];
  int *hadapt = (int *)data[1];

  user_data_t *qdata = info->quad->p.user_data;
  qdata->iref = hadapt[*k];

  ++(*k);
}

void P2N_FUN(mark_elements)(int *hadapt)
{
  p4est_t *p4est = stored_p4est;
  int k = 0;
  void *data[] = {&k, hadapt};
  p4est_iterate(p4est, NULL, data, mark_element_fun, NULL, NULL
#ifdef P4_TO_P8
                ,
                NULL
#endif
                );
}

#define P2N_FLAG_COARSEN (1 << 0)
#define P2N_FLAG_REFINE (1 << 1)

static void quadrant_replace(p4est_t *p4est, p4est_topidx_t which_tree,
                             int num_outgoing, p4est_quadrant_t *outgoing[],
                             int num_incoming, p4est_quadrant_t *incoming[])
{
  if (num_outgoing == 1)
  {
    /* Refining: copy data to all children */
    for (int c = 0; c < num_incoming; ++c)
    {
      user_data_t *in_ud = incoming[c]->p.user_data;
      user_data_t *out_ud = outgoing[0]->p.user_data;
      in_ud->iref = out_ud->iref;
    }
  }
  else
  {
    user_data_t *in_ud = incoming[0]->p.user_data;
    user_data_t *out_ud = outgoing[0]->p.user_data;
    in_ud->iref = out_ud->iref;
  }
}

static void adapt_p4est(uint8_t flags, int32_t *num_dst)
{
  p4est_t *p4est = stored_p4est;

  if (flags & P2N_FLAG_COARSEN)
    p4est_coarsen_ext(p4est, 0, 0, coarsen_list_fn, NULL, quadrant_replace);

  if (flags & P2N_FLAG_REFINE)
    p4est_refine_ext(p4est, 0, -1, refine_list_fn, NULL, quadrant_replace);

  p4est_balance_ext(p4est, P4EST_CONNECT_FULL, NULL, quadrant_replace);

  *num_dst = p4est->local_num_quadrants;
}

void P2N_FUN(coarsen_p4est)(int32_t *num_dst)
{
  adapt_p4est(P2N_FLAG_COARSEN, num_dst);
}

void P2N_FUN(refine_p4est)(int32_t *num_dst)
{
  adapt_p4est(P2N_FLAG_REFINE, num_dst);
}

void P2N_FUN(coarsen_refine_p4est)(int32_t *num_dst)
{
  adapt_p4est(P2N_FLAG_COARSEN | P2N_FLAG_REFINE, num_dst);
}

/* Sohail Reddy (08/05/2020): comment out to use with new p4est
void P2N_FUN(dump_forest)(const char *fname)
{
  p4est_vtk_write_all(stored_p4est, NULL, 1, 1, 1, 1, 0, 0, 0, fname);
}
*/

static void get_elm_lvl(p4est_iter_volume_info_t *info, void *user_data)
{
  void **data = user_data;

  int *k = (int *)data[0];
  int8_t *lvl = (int8_t *)data[1];

  lvl[*k] = info->quad->level;

  ++(*k);
}
void P2N_FUN(get_element_lvl)(int8_t *lvl, int32_t *N)
{
  p4est_t *p4est = stored_p4est;

  if (*N != p4est->local_num_quadrants)
    SC_ABORT("invalid array size");

  int k = 0;
  void *data[] = {&k, lvl};
  p4est_iterate(p4est, NULL, data, get_elm_lvl, NULL, NULL
#ifdef P4_TO_P8
                ,
                NULL
#endif
                );
}

void P2N_FUN(repartition)(int64_t *qid_src, int64_t *qid_dst)
{
  p4est_t *p4est = stored_p4est;
  for (int k = 0; k < p4est->mpisize + 1; ++k)
    qid_src[k] = p4est->global_first_quadrant[k];

  p4est_partition(p4est, 1, NULL);

  for (int k = 0; k < p4est->mpisize + 1; ++k)
    qid_dst[k] = p4est->global_first_quadrant[k];
}

/* END OF BFAM CODE */



/*
 * Read GMSH grids
 */
int READ_MESH(void)
{
   
#define P4_TO_P8
	    
    int *is_cube=1;
    int *nnx=10;
    int *nny=10;
    int *nnz=10;
    int *nref_levs=0;
    int *read_external_grid_flg=1;
    int *xperiodic_flg=0;
    int *yperiodic_flg=0;
    int *zperiodic_flg=0;
#ifdef P4_TO_P8
    printf(" DEFINED XXXXXXx\n");

    //start p4est
    int p4est_log_level_c;
    p4est_log_level_c = 6;
    p4est_start_(p4est_log_level_c);
#endif
    /*p8esttonuma_init(is_cube, nnx, nny, nnz, nref_levs, read_external_grid_flg, \
      xperiodic_flg, yperiodic_flg, zperiodic_flg);
    */
    /*    p8esttonuma_fill_data(nop, xgl, is_dg, iboundary, read_external_grid_flg, p2n, \
	  lrestoring_sponge);
    */
    /* p8esttonuma_get_mesh_scalars(p2n, npoin_cg, nelem, num_nbh, \
       num_send_recv_total, nbsido, nface, nboun, nNC)*/
	    
    /*
    
    int *refine_coarsen_elements;
    refine_coarsen_elements = (int *)malloc(1*sizeof(int *));
    int read_external_grid_flg = 1;
    int xperiodic_flg = 0;
    int yperiodic_flg = 0;
    int zperiodic_flg = 0;
    int is_cube;
    int nnx, nny, nnz;
    int orient = 0;
    int nref_levs = 0;
    int refinement_levels_h = 0;

    printf(" XXX\n");
    is_cube = 1;
    nnx = 19;
    nny = 19;
    nnz = 19;
    p8esttonuma_init(is_cube, nnx, nny, nnz, nref_levs,		\
		     read_external_grid_flg,			\
		     xperiodic_flg, yperiodic_flg, zperiodic_flg);

    printf(" YYY\n");
    /* ! TODO: fix to use xglx, xgly, xglz
       call p8esttonuma_fill_data(nop, xgl, is_dg, iboundary, read_external_grid_flg, p2n, lrestoring_sponge)

       call p8esttonuma_get_mesh_scalars(p2n, npoin_cg, nelem, num_nbh, &
       num_send_recv_total, nbsido, nface, nboun, nNC)
    */
    /*
        free(refine_coarsen_elements);
    */
    return 0;
}
