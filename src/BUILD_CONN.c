/*
 * This function builds the connectivity matrix:
 *
 * Simone Marras, Aug 2012
 */
#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

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
