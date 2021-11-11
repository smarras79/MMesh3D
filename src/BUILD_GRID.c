/*
 * This function computes the grid coordinates according to the type defined by the user in input:
 *
 * Simone Marras, Aug 2012
 *
 */
#include <stdlib.h>
#include <string.h>

#include "BUILD_GRID.h"
#include "GLOBAL_VARS.h"
#include "GRID_COORD.h"
#include "MEMORY.h"
#include "MYDEFINE.h"
#include "SURFACES.h"

int BUILD_GRID(void)
{
  if( !strcmp(problem[0], "rectangle") || !strcmp(problem[0], "cube") )

    if( nop == 1)
      GRID_COORD();
    else
      {  
	printf("\n");
	printf(" ERROR: HIGH-ORDER GRID: NOT PROGRAMMED YET!!!\n"); 
	printf("        The program will EXIT now\n");
	printf("\n");
	exit(1);
	
	/*SURFACES_HIGH_ORDER(1, nnodesx, nnodesy, nnodesz, nop, ngl, nelem, \
	 *		  BDY_COORDS, problem, parameters, COORDS, rarray, deltaLon, deltaLat);
	 */
      }
  else if( !strncmp(problem[0], "topograpny", 4) || \
	   !strncmp(problem[0], "TOPOGRAPHY", 4) || \
	   !strncmp(problem[0], "mountain", 4)   || \
	   !strncmp(problem[0], "agnesi", 4))
    {

      if( nop == 1)
	{
	  SURFACES(1, nnodesx, nnodesy, nnodesz, nelem, BDY_COORDS, problem, parameters, 
		   COORDS, rarray, deltaLon, deltaLat);
	}
      else
	{
	  printf("\n");
	  printf(" ERROR: HIGH-ORDER GRID: NOT PROGRAMMED YET!!!\n"); 
	  printf("        The program will EXIT now\n");
	  printf("\n");
	  exit(1);
	  /*  SURFACES_HIGH_ORDER(1, nnodesx, nnodesy, nnodesz, nop, ngl, nelem, BDY_COORDS, problem, parameters, COORDS, rarray, deltaLon, deltaLat);
	   */
	}
      
      /*
       * Free the memory for 'rarray' because it won't be used anymore
       */
      if( !strcmp(problem[9], "NCDF")    ||		\
	  !strcmp(problem[9], "ncdf")    ||		\
	  !strcmp(problem[9], "netcdf")  ||		\
	  !strcmp(problem[9], "NetCDF")  ||		\
	  !strcmp(problem[9], "Ncdf"))
	{
	  printf(" #------------------------------------------------------------------#\n");
	  printf(" #------------------------------------------------------------------#\n");
	  //free_dvector(rarray, 0,nnodesx*nnodesy);
	}
      else if( !strcmp(problem[9], "text")   ||	\
	       !strcmp(problem[9], "txt")    ||	\
	       !strcmp(problem[9], "noaa")   ||	\
	       !strcmp(problem[9], "NOAA")   || \
	       !strcmp(problem[9], "DEM")    || \
	       !strcmp(problem[9], "Dem")    || \
	       !strcmp(problem[9], "dem"))
	{
	  MEMORY_DEALLOCATE(3);
	}
    }
  
  return 0;
}
