/*******************************************************************************
 *
 * READTOPOtxt.c
 *
 * This function reads the topography form a file of shape [1:nnodes][3]
 * where the first and second column are the ordered lat-lon coordinates
 * and the third column is the height of topography at that specific
 * coordinate point.
 *
 * 1) XYZ files from NOAA
 *  
 * READTOPOtxt_header() reads the parameters from the header file (*.hdr)
 * READTOPOtxt_file()   reads the actual file of coordinates      (*.xyz)
 *
 * 2) DEM files
 * 
 * READTOPO_DEM()       reads a DEM file from NOAA page (file extension: *.asc)
 *
 ********************************************************************************/
#include "myinclude.h"

int READTOPOtxt_header(char *txt_inputfile, 
		       int *nlon, int *nlat,
		       double *deltaLon, double *deltaLat)
{

  FILE *file_ID;
  char *p;
  
  double lonmin,lonmax,latmin,latmax,deltacell;
  
  char word[64];  

  //Add the extension of the file name:
  txt_inputfile = strtok (txt_inputfile,".");
  strcat(txt_inputfile, ".hdr");

  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  
  if((file_ID = fopen(txt_inputfile, "r")) == NULL)
    printf("[READTOPOtxt_header] The file %s could not be open\n", txt_inputfile);
  else{
    printf(" # Opening header file %s ........... DONE\n", txt_inputfile);

    fscanf(file_ID, "%s %d\n", word, nlon);        //NCOLS (LON)
    fscanf(file_ID, "%s %d\n", word, nlat);        //NROWS (LAT)
    fscanf(file_ID, "%s %lf\n", word, &lonmin);    //XLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &latmin);    //YLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &deltacell); //CELLSIZE
    
    lonmax = lonmin + *nlon*deltacell;
    latmax = latmin + *nlat*deltacell;
    
  }
  
  *deltaLon = lonmax - lonmin;
  *deltaLat = latmax - latmin;

  return 0;
}

int READTOPOtxt_file(char *txt_inputfile, char *topoBathy_flg, 
		     double *rarray, int nnodesx, int nnodesy)
{
  
  FILE *file_ID;
  
  int i,j,k;
  int nnodes;
  char word[64];

  //  double LonLatZ[nnodesx][nnodesy], dummy;
  double dummy;
  double LonLatZ[nnodesx][nnodesy];

  //Add the extension of the file name:
  txt_inputfile = strtok (txt_inputfile,".");
  strcat(txt_inputfile, ".xyz");

  if((file_ID = fopen(txt_inputfile, "r")) == NULL)
    printf("The file %s could not be open\n", txt_inputfile);
  else{
    k=0;
    while(!feof(file_ID)){
      k = k + 1;
      fscanf(file_ID, "%lf %lf %lf\n", &dummy, &dummy, &rarray[k]);
    }
  } 

  k = 0;
  for(j=nnodesy-1; j>=0; j--){
    for(i=0; i<nnodesx; i++){
      k = k + 1;
      LonLatZ[i][j] = rarray[k];
      //printf(" k = %d lollatz = %lf\n", k,LonLatZ[i][j]);
    }
  }
  
  k = 0;

  for(j=0; j<nnodesy; j++){    
    for(i=0; i<nnodesx; i++){
   
      k = k + 1;
      
      if( !strncmp(topoBathy_flg, "topography", 4))
	{
	  /*
	   * the keyword in input "topo" means that only the positive Z are read
	   */
	  if(  LonLatZ[i][j] < 0.0 )
	    rarray[k] = 0.0;
	  else
	    rarray[k] = LonLatZ[i][j];
	}
      else if( !strncmp(topoBathy_flg, "bathymetry", 4))
	{
	  /*
	   * the keyword in input "bathy" means that only the positive Z are read
	   */
	  if(  LonLatZ[i][j] > 0.0 )
	    rarray[k] = 0.0;
	  else
	    rarray[k] = LonLatZ[i][j];
	}
      else if( !strcmp(topoBathy_flg, "all"))
	{
	  /*
	   * the keyword in input "bathy" means that only the positive Z are read
	   */
	  rarray[k] = LonLatZ[i][j];
	}
      else
	{
	  printf("\n # ERROR MESSAGE\n");
	  printf(" # [READTOPOtxt_file] ERROR IN INPUT:\n");
	  printf(" #                    use one of:\n");
	  printf(" #                    - topography\n");
	  printf(" #                    - bathymetry\n");
	  printf(" #                    - all\n");
	  printf(" # The program will EXIT now !!! \n");
	  printf(" # Open the Input file and fix it and run the code again\n\n");
	  exit(1);
	}
      
    }
  }

  return 0;
}

/*********************************************************************************
 *
 * deg2meters.c
 *
 * This function reads in a differential dimension in 
 * degrees (deltalat and deltalon), and return the equivalent in meters
 *
 * We use the approx.: 1deg <--> 111.3km (ONLY VALID AT THE EQUATOR THEN!)
 * 
 *********************************************************************************/
int deg2meters(double deltaLon, double deltaLat, double *deltaX, double *deltaY)
{
  
  *deltaX = deltaLon*111.3*1000.0;
  *deltaY = deltaLat*111.3*1000.0;

  return 0;
}


/**********************************************************************************
 *
 * READTOPO_DEM.c
 *
 * This function reads a DEM file (*.asc) as it comes from the NOAA
 * repository (http://www.ngdc.noaa.gov/dem/selectdem.jsp?regionName=West%20Coast)
 * 
 **********************************************************************************/
int READTOPO_DEM_header(char *txt_inputfile, 
			int *nlon, int *nlat,
			double *deltaLon, double *deltaLat)
{

  FILE *file_ID;
  char *p;
  
  double lonmin,lonmax,latmin,latmax,deltacell;
  int nodata_value,k;

  char word[64];  

  //Add the extension of the file name:
  txt_inputfile = strtok (txt_inputfile,".");
  strcat(txt_inputfile, ".asc");

  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  
  if((file_ID = fopen(txt_inputfile, "r")) == NULL)
    printf("[READTOPO_DEM_header] The file %s could not be open\n", txt_inputfile);
  else{
    printf(" # Opening DEM file %s ........... DONE\n", txt_inputfile);

    fscanf(file_ID, "%s %d\n", word, nlon); //NCOLS (LON)
    fscanf(file_ID, "%s %d\n", word, nlat); //NROWS (LAT)
    fscanf(file_ID, "%s %lf\n", word, &lonmin); //XLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &latmin); //YLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &deltacell); //CELLSIZE
    
    lonmax = lonmin + *nlon*deltacell;
    latmax = latmin + *nlat*deltacell;
    
  }
  
  *deltaLon = lonmax - lonmin;
  *deltaLat = latmax - latmin;

  return 0;
}


int READTOPO_DEM_file(char *txt_inputfile, char *topoBathy_flg, 
		      double *rarray, int nlon, int nlat)
{

  FILE *file_ID;
  char *p;  

  double lonmin,lonmax,latmin,latmax,deltacell,dummy;
  int nodata_value,k;

  char word[64];

  //Add the extension of the file name:
  txt_inputfile = strtok (txt_inputfile,".");
  strcat(txt_inputfile, ".asc");

  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  
  if((file_ID = fopen(txt_inputfile, "r")) == NULL)
    printf("[READTOPO_DEM_header] The file %s could not be open\n", txt_inputfile);
  else{
    printf(" # Opening DEM file %s ........... DONE\n", txt_inputfile);

    fscanf(file_ID, "%s %d\n", word, &nlon);          //NCOLS (LON)
    fscanf(file_ID, "%s %d\n", word, &nlat);          //NROWS (LAT)
    fscanf(file_ID, "%s %lf\n", word, &lonmin);       //XLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &latmin);       //YLLCENTER
    fscanf(file_ID, "%s %lf\n", word, &deltacell);    //CELLSIZE
    fscanf(file_ID, "%s %d\n", word, &nodata_value); //NODATA_value
    
    k=0;
    while(!feof(file_ID)){
      k = k + 1;
      fscanf(file_ID, "%lf \n", &dummy);
      rarray[k] = dummy;
    }

    printf(" # Reading DEM file %s ........... DONE\n", txt_inputfile);
  }

  return 0;
}
