/*
 * This functions reads the topography files (type of topography defined in input by the user):
 * 
 * Simone Marras, Aug 2012
 */
#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

int READ_TOPOGRAPHY(void)
{

  // The file to be read is a Netcdf:
  if( !strncmp(problem[0], "topograpny", 4) || !strncmp(problem[0], "TOPOGRAPHY", 4) || \
      !strncmp(problem[0], "Topography", 4) || !strncmp(problem[0], "mountain", 4)   || \
      !strncmp(problem[0], "orography", 4) || !strncmp(problem[0], "OROGRAPHY", 4) )
    {
      if( !strcmp(problem[9], "NCDF")   ||	\
	  !strcmp(problem[9], "ncdf")   ||	\
	  !strcmp(problem[9], "netcdf") ||	\
	  !strcmp(problem[9], "NetCDF") ||	\
	  !strcmp(problem[9], "Ncdf"))
	{
	  printf(" #------------------------------------------------------------------#\n");
	  printf("   !!! NETCDF REDING NOT PROGRAMMED YET\n");
	  printf("   !!! ERROR in %s: TOPOGRAPHY must be either noaa or FUNCTION !!! \n", inputfile);
	  printf("   !!! open %s and change it with noaa (if you have ancdf file) or FUNCTION\n", problem[9]);
	  printf("   !!! The program will exit now\n");
	  printf(" #------------------------------------------------------------------#\n");
	}
      else if( !strcmp(problem[9], "DEM")     ||	\
	       !strcmp(problem[9], "Dem")     ||	\
	       !strcmp(problem[9], "dem"))
	{
	
	  /* Get the size of the topography (NOAA) file to be read */
	  READTOPO_DEM_header(problem[5], &nnodesx, &nnodesy, &deltaLon, &deltaLat);

	  /* GEOGRID rarry */
	  nnodes = nnodesx*nnodesy;
	  MEMORY_ALLOCATE(3);
	  
	  READTOPO_DEM_file(problem[5], problem[11], rarray, nnodesx, nnodesy);

	}
      else if( !strcmp(problem[9], "text")   ||	\
	       !strcmp(problem[9], "txt")    ||	\
	       !strcmp(problem[9], "noaa")   ||	\
	       !strcmp(problem[9], "NOAA"))
	{

	  /* Get the size of the topography (NOAA) file to be read */
	  READTOPOtxt_header(problem[5], &nnodesx, &nnodesy, &deltaLon, &deltaLat);

	  /* GEOGRID rarry */
	  nnodes = nnodesx*nnodesy;
	  MEMORY_ALLOCATE(3);
	
	  READTOPOtxt_file(problem[5], problem[11], rarray, nnodesx, nnodesy);
	}
    }

  //If the topography is read from file, adjust the number of elements and nodes
  //according to the values of nnodesx,y,z read from the topography file:
  nnodes = nnodesx*nnodesy*nnodesz;
  nelx = nnodesx-1;
  nely = nnodesy-1;
  nelz = nnodesz-1;

  //Number of nodes for high-order elements (nop>1):
  ngl=nop + 1;
  nx=nelx*nop + 1;
  ny=nely*nop + 1;
  nz=nelz*nop + 1;
    
  //Global sizes
  nnodes= nx*ny*nz;
  nelem = nelx*nely*nelz;
  nboun = 2*nelx*nely + 2*nelx*nelz + 2*nely*nelz;
  
  if( !strncmp(problem[2],"WEDGE",4) )
    nelem = 2*nelem;


  return 0;
}
