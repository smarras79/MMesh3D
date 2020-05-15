
#include "myinclude.h"

#define RANK_time 1
#define MAX_NVARS 100
#define VARNAME_LENGHT 24
#define DUMB 10
#define TIMES 0
#define LONS  0
#define LATS  0

#define flg_dummy 2

int GETNCDF_SIZE(char *ncdf_inputfile, int *nlon, int *nlat)
{
  FILE *mat_id;
  
  const char *ncdf_filename[128];
  
  int i,j,k,l, dumb;
  
  int status;			//Error status
  int ncid;			//NetCDF ID
  int time_id;			//Variable ID
 
  int ndims, nvars, ngatts, unlimdimid;	//Existing variables characteristics
  
  static char varname[NC_MAX_NAME+1];
  char unlmtd_name[NC_MAX_NAME+1];
   
  int  n;
  int latid, lonid, zid, recid;
  size_t t, time, nvlevels, unlimtd_length;
  size_t lonlength, latlength;

  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  *ncdf_filename = ncdf_inputfile;
  
  //Take the terminal input of the path and file name
  printf("\n - f_name = %s\n", *ncdf_filename);
 
  //Open dataset to access it in read-only (putting 0 instead of NC_WRITE) or read-write mode:
  status = nc_open(*ncdf_filename, NC_NOWRITE, &ncid);
 
  status = nc_inq(ncid, &ndims,  &nvars, &ngatts, &unlimdimid);
  printf("\n - File: %s\n\n - ID = %d\n - ndims = %d\n - nvars = %d\n - ngatts = %d\n\n", \
	 *ncdf_filename, ncid, ndims, nvars, ngatts);
  
  //Read values for "dimensions":
  
  // get ID and NAME (unlmtd_name) of unlimited dimension
  status = nc_inq_unlimdim(ncid, &recid);
  status = nc_inq_dim(ncid, recid, unlmtd_name, &unlimtd_length);
  printf(" - Unlimited Dim. Name : %s \n - Unlimited Dim. Length : %d\n", unlmtd_name, unlimtd_length);
 

  /* Set #define flg_dummy = 1 for USA file
   *                       = 2 for SPAIN file
   *                       = 3 for WORLD file
   * for testing only.
   */
  
  switch (flg_dummy)
    {
    case 1:
      //NCsouth-eastUSA.nc
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "west_east", &lonid);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "south_north", &latid);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      printf(" - south_north (LAT) : %d\n", latlength);
      
      //n grid points in the vertical (z) direction
      status = nc_inq_dimid(ncid, "bottom_top", &zid);
      status = nc_inq_dimlen(ncid, zid, &nvlevels);
      printf(" - bottom_top (nvlevels) : %d\n", nvlevels);
      break;
      
    case 2:
      //SPAIN
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "ROW", &lonid);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "COL", &latid);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      printf(" - south_north (LAT) : %d\n", latlength);
      
      //n grid points in the vertical (z) direction
      status = nc_inq_dimid(ncid, "LAY", &zid);
      status = nc_inq_dimlen(ncid, zid, &nvlevels);
      printf(" - bottom_top (nvlevels) : %d\n", nvlevels);
      break;
      
    case 3:
      //WORLD
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "Lon-000", &lonid);
      //handle_error(status);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      //handle_error(status);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "Lat-000", &latid);
      //handle_error(status);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      //handle_error(status);
      printf(" - south_north (LAT) : %d\n", latlength);
      break;
    }


    
  //Dynamic allocation:
  //Only take the levels that are physically existent, and not only numerically:
  //this means, for example, that if latlength = 2 it means that the dimension is
  //only numerical in the sense that WRF needs 2 grid points to compute the flux.
  //In this case, we don't consider that value to compute the total size of the problem.
  
  if( lonlength < 4)
    lonlength = 1;  
  if( latlength < 4)
    latlength = 1;
  
  //FOR NOW WE CARE ONLY ABOUT THE FIRST LEVEL and don't care about any other:
  //if( nvlevels < 3 ){
    nvlevels = 1;
    //  printf(" - NVLEVELS %d\n", nvlevels);
    //}
  
  *nlon = lonlength;
  *nlat = latlength;

  return 0;
}

int READNCDF(char *ncdf_inputfile, float *rarray)
{
  FILE *mat_id;
  
  const char *ncdf_filename[128];
  
  int i,j,k,l, dumb;
  
  int status;			//Error status
  int ncid;			//NetCDF ID
  int time_dim;                 //Dimension ID
  int time_id;			//Variable ID
  int time_dimids[0];		//Variable shape
  int ndims, nvars, ngatts, unlimdimid;	//Existing variables characteristics
  size_t time_len = NC_UNLIMITED;
  
  int n;
  int rh_id;                         /* variable ID */
  int rh_ndims, ndimsp;               /* number of dims */
  int var_id;
  int latid, lonid, zid, recid;
  size_t t, time, latlength, lonlength, nvlevels, unlimtd_length;
  
  int vardim[] = {1};           //Now only for Temperature (scalar)
  int *dims;
  int *vardims;

  //For node centered variables (it'll be 1 for all! Never cell centered for our purposes)
  int centering[] = {1}; 
  
  static char varname[NC_MAX_NAME+1];
  char unlmtd_name[NC_MAX_NAME+1];
  
  char **varnames;
  char **varnames2;
  
  float **vars;
  float *var_vals;
  
  int npts, ipoin;
  int ncX, ncY, ncZ;
  int ncells;
  
  //Attributes:
  nc_type vr_type, t_type;   /* attribute types */
  size_t  vr_len, t_len;     /* attribute lengths */
  
  /* Rectilinear mesh coordinates. */
  size_t size;
  
  float DX, DY, z_scale;
  int KDE;		 //VETICAL LEVELS. kde is used in WRF so we keep the same notation.
  // This value is given in the namelist.input as: "e_vert"
  

  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  *ncdf_filename = ncdf_inputfile;
  
  //Take the terminal input of the path and file name
  printf("\n - File_name = %s\n", *ncdf_filename);
 
  //Open dataset to access it in read-only 
  //(putting 0 instead of NC_WRITE) or read-write mode:
  status = nc_open(*ncdf_filename, NC_NOWRITE, &ncid);
  
  status = nc_inq(ncid, &ndims,  &nvars, &ngatts, &unlimdimid);
  printf("\n - File: %s\n\n - ID = %d\n - ndims = %d\n - nvars = %d\n - ngatts = %d\n\n", \
	 *ncdf_filename, ncid, ndims, nvars, ngatts);
  
  // get ID and NAME (unlmtd_name) of unlimited dimension
  status = nc_inq_unlimdim(ncid, &recid);
  status = nc_inq_dim(ncid, recid, unlmtd_name, &unlimtd_length);
  printf(" - Unlimited Dim. Name : %s \n - Unlimited Dim. Length : %d\n", unlmtd_name, unlimtd_length);
  
  /* Set #define flg_dummy = 1 for USA file
   *                       = 2 for SPAIN file
   *                       = 3 for WORLD file
   * for testing only.
   */
  
  switch (flg_dummy)
    {
    case 1:
      //NCsouth-eastUSA.nc
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "west_east", &lonid);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "south_north", &latid);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      printf(" - south_north (LAT) : %d\n", latlength);
      
      //n grid points in the vertical (z) direction
      status = nc_inq_dimid(ncid, "bottom_top", &zid);
      status = nc_inq_dimlen(ncid, zid, &nvlevels);
      printf(" - bottom_top (nvlevels) : %d\n", nvlevels);
      break;
      
    case 2:
      //SPAIN
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "ROW", &lonid);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "COL", &latid);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      printf(" - south_north (LAT) : %d\n", latlength);
      
      //n grid points in the vertical (z) direction
      status = nc_inq_dimid(ncid, "LAY", &zid);
      status = nc_inq_dimlen(ncid, zid, &nvlevels);
      printf(" - bottom_top (nvlevels) : %d\n", nvlevels);
      break;
      
    case 3:
      //WORLD
      //n grid points in the long. direction
      status = nc_inq_dimid(ncid, "Lon-000", &lonid);
      //handle_error(status);
      status = nc_inq_dimlen(ncid, lonid, &lonlength);
      //handle_error(status);
      printf(" - west_east (LON) : %d\n", lonlength);
      
      //n grid points in the lat. direction
      status = nc_inq_dimid(ncid, "Lat-000", &latid);
      //handle_error(status);
      status = nc_inq_dimlen(ncid, latid, &latlength);
      //handle_error(status);
      printf(" - south_north (LAT) : %d\n", latlength);
      break;
    }
 
 
  //Dynamic allocation:
  //Only take the levels that are physically existent, and not only numerically:
  //this means, for example, that if latlength = 2 it means that the dimension is
  //only numerical in the sense that WRF needs 2 grid points to compute the flux.
  //In this case, we don't consider that value to compute the total size of the problem.
  
  if( lonlength < 4)
    lonlength = 1;  
  if( latlength < 4)
    latlength = 1;
  
  //FOR NOW WE CARE ONLY ABOUT THE FIRST LEVEL and don't care about any other:
  //if( nvlevels < 3 || nvlevels >= 10E+03){
    nvlevels = 1;
    //  printf(" - NVLEVELS %d\n", nvlevels);
    //}

  // Dynamic allocation of the arrays to pointers char **varnames and float **vars;
  // This line only allocates the number of entries but NOT the length of the entries:
  dims = (int*) malloc(3*sizeof(int *));
  vardims = (int*) malloc(nvars* sizeof(int *));

  varnames = (char**) malloc(nvars* sizeof(char *));
  for(i=0; i<nvars; i++)
    varnames[i] = (char*) malloc( VARNAME_LENGHT * sizeof(char *));
  
  varnames2 = (char**) malloc( 1 * sizeof(char *));
  varnames2[0] = (char*) malloc( VARNAME_LENGHT * sizeof(char *));
  
  vars = (float**) malloc(nvars* sizeof(float *));
  for(i=0; i<nvars; i++)
    vars[i] = (float*) malloc(DUMB * sizeof(float *));

  //Fill in the dims array: dims[] = {NX, NY, NZ}
  dims[0] = lonlength;
  dims[1] = latlength;
  dims[2] = nvlevels;
  size = lonlength*latlength*nvlevels;

  var_vals = (float *) malloc(size* sizeof(float *));

  //Allocate and read the ncdf variable values:
  n = 0;
  while(n < nvars)
    {
      
      vardims[n] = 1; //This must be changed as to be assigned if vector (3) or scalar (1)
      var_id = n;
      
      // Allocation of the length of each entry:
      status = nc_inq_varname(ncid, n, varname);
      printf("Variable: %s\n", varname);
      
      varnames[n] = (char*) malloc( sizeof(varname) * sizeof(char *));
      
      //Assign each entry of varnames:
      strcpy(varnames[n], varname);
            
      if( !strcmp(varnames[n],"PB") ){
	strcpy(varnames2[0], varnames[n]);
	printf(" Variable name: %s\n", varnames2[0]);
	var_id = n;
	
	nc_inq_varndims (ncid, var_id, &ndimsp);
	printf("Varname[%d]: %s -> var_id: %d -> dim : %d\n", n, varnames[n], var_id, ndimsp);	
	status = nc_get_var_float(ncid, var_id, var_vals);
      }
      else if( !strcmp(varnames[n],"HGT_M") || !strcmp(varnames[n],"HT") ){
	strcpy(varnames2[0], varnames[n]);
	printf(" Variable name: %s\n", varnames2[0]);
	var_id = n;
	
	nc_inq_varndims (ncid, var_id, &ndimsp);
	printf("Varname[%d]: %s -> var_id: %d -> dim : %d\n", n, varnames[n], var_id, ndimsp);	
	status = nc_get_var_float(ncid, var_id, var_vals);
	
      }
      else if( !strcmp(varnames[n],"GMAO-2D__PHIS") ){
	strcpy(varnames2[0], varnames[n]);
	printf(" Variable name: %s\n", varnames2[0]);
	var_id = n;
	
	nc_inq_varndims (ncid, var_id, &ndimsp);
	printf("Varname[%d]: %s -> var_id: %d -> dim : %d\n", n, varnames[n], var_id, ndimsp);	
	status = nc_get_var_float(ncid, var_id, var_vals);
      }

      n = n + 1;
    }

  ipoin=0;
  for(j=0; j<lonlength; j++){
    for(k=0; k<latlength; k++){
      rarray[ipoin+1] =  var_vals[ipoin];
      //printf("level=%d k=%d  rhvals[%d]: %f %f\n",i, ipoin, ipoin, var_vals[ipoin],rarray[ipoin+1] );
      ipoin++;
    }
  }
    
  
  //Close the file and write it:
  status = nc_close(ncid);
  
    //Free memory:
  free(dims);
  free(vardims);
  free(var_vals);
  for(i=0; i<nvars; i++){
    free(varnames[i]);
  }
  free(varnames);
  free(varnames2[0]);
  free(varnames2);
  for(i=0; i<nvars; i++)
    free(vars[i]);
  free(vars);	
  
  return 0;
}


/*Code to open, read and write into a txt file a matrix variable from a NetCDF file.
 *This routine was specifically written in the context of topography reading.
 *
 * Author: Simone Marras
 * Created on: 24 May 2008
 * simone.marras@gmail.com
 *
 *************************************************************************************/

int TOPOfromNCDF(char *ncdf_filename, char *input_var, double *rarray)
{
	
  int i,j,k;
		
  unsigned int wrt2file = 0;
  int ncid;     	                //NetCDF ID
  int input_var_id; 	                //variable ID
	 
  size_t s1_cnt, s2_cnt;
	
  int status;				//Error status
	
  int time_dim;				//Dimension ID
  int time_dims[RANK_time];
  int time_id;				//Variable ID
  int time_dimids[0];			//Variable shape
  int ndims, nvars, ngatts, unlimdimid;	//Existing variables characteristics
  size_t time_len = NC_UNLIMITED;
	
  static size_t start[] = {0, 0, 0}; /* start at first value */
  static size_t count[] = {TIMES, LONS, LATS};
  int VAR_INT[LONS*LATS]; /* array to hold values */
  double VAR_DBL[LONS*LATS];

  FILE *mat_id;
	
  /***********************************************************
   * Preliminary computations to prepare some strings
   ***********************************************************/
  
  //Open VAR_INTset to access it in read-only (putting 0 instead of NC_WRITE) or read-write mode:
  status = nc_open(ncdf_filename, NC_NOWRITE, &ncid);
  if(status != NC_NOERR){
    printf("\n\t!ERROR in Opening file '%s'\n\t!The program will now exit. \n\t!Check that the file you want to open exists in your directory!\n\n", ncdf_filename);
    return 1;}
  else printf("\nOK - File %s opened with ID: %d\n", ncdf_filename, ncid);
  
  status = nc_inq(ncid, &ndims,  &nvars, &ngatts, &unlimdimid);
  
  status = nc_inq_varid (ncid, input_var, &input_var_id);
  if (status != NC_NOERR){
    printf("\nError: you asked for a variable that is NOT is you *.nc file\nCheck and run the program again\n\n");
    return 1;
  }
  
  /* read values from netCDF variable */
  if(strcmp(input_var,"LANDMASK") == 0){
    status = nc_get_vara_int(ncid, input_var_id, start, count, VAR_INT);
    //printf("status = %d\n", status);
  }
  else if(strcmp(input_var,"SLPX") == 0 || strcmp(input_var,"HGT_M") == 0){
    status = nc_get_vara_double(ncid, input_var_id, start, count, rarray);
    //printf("status = %d\n", status);
  }
  
  return 0;
}


/*********************************************************************
 *
 * READTOPOtxt.c
 *
 * This function reads the topography form a file of shape [1:nnodes][3]
 * where the first and second column are the ordered lat-lon coordinates
 * and the third column is the height of topography at that specific
 * coordinate point.
 * 
 * READTOPOtxt_header() reads the parameters from the header file
 * READTOPOtxt_file()   reads the actual file of coordinates
 *
 *
 *********************************************************************/
int READTOPOtxt_header(char *txt_inputfile, 
		       int *nlon, int *nlat,
		       float *deltaLon, float *deltaLat)
{

  FILE *file_ID;
  char *p;
  
  float lonmin,lonmax,latmin,latmax,deltacell;
  
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

    fscanf(file_ID, "%s %d\n", word, nlon); //NCOLS (LON)
    fscanf(file_ID, "%s %d\n", word, nlat); //NROWS (LAT)
    fscanf(file_ID, "%s %f\n", word, &lonmin); //XLLCENTER
    fscanf(file_ID, "%s %f\n", word, &latmin); //YLLCENTER
    fscanf(file_ID, "%s %f\n", word, &deltacell); //CELLSIZE
    
    lonmax = lonmin + *nlon*deltacell;
    latmax = latmin + *nlat*deltacell;
    
  }
  
  *deltaLon = lonmax - lonmin;
  *deltaLat = latmax - latmin;

  return 0;
}

int READTOPOtxt_file(char *txt_inputfile, char *topoBathy_flg, 
		     float *rarray, int nnodesx, int nnodesy)
{
  
  FILE *file_ID;
  
  int i,j,k;
  int nnodes;
  char word[64];

  //  float LonLatZ[nnodesx][nnodesy], dummy;
  float dummy;
  float LonLatZ[nnodesx][nnodesy];

  //Add the extension of the file name:
   //Add the extension of the file name:
  txt_inputfile = strtok (txt_inputfile,".");
  strcat(txt_inputfile, ".xyz");

  if((file_ID = fopen(txt_inputfile, "r")) == NULL)
    printf("The file %s could not be open\n", txt_inputfile);
  else{
    k=0;
    while(!feof(file_ID)){
      k = k + 1;
      fscanf(file_ID, "%f %f %f\n", &dummy, &dummy, &rarray[k]);
    }
  } 

  k = 0;
  for(j=nnodesy-1; j>=0; j--){
    for(i=0; i<nnodesx; i++){
      k = k + 1;
      LonLatZ[i][j] = rarray[k];
      //printf(" k = %d lollatz = %f\n", k,LonLatZ[i][j]);
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

/*********************************************************************
 *
 * deg2meters.c
 *
 * This function reads in a differential dimension in 
 * degrees (deltalat and deltalon), and return the equivalent in meters
 *
 * We use the approx.: 1deg <--> 111km
 * 
 *********************************************************************/
int deg2meters(float deltaLon, float deltaLat, float *deltaX, float *deltaY)
{
  
  *deltaX = deltaLon*111*1000.0;
  *deltaY = deltaLat*111*1000.0;
 
  return 0;
}
