#ifndef TOPOFROMTXT_H
#define TOPOFROMTXT_H

int READTOPOtxt_header(char *txt_inputfile, 
		       int *nlon, int *nlat,
		       double *deltaLon, double *deltaLat);

int READTOPOtxt_file(char *txt_inputfile, char *topoBathy_flg, 
		     double *rarray, int nnodesx, int nnodesy);

int deg2meters(double deltaLon, double deltaLat, double *deltaX, double *deltaY);

int READTOPO_DEM_header(char *txt_inputfile, 
		       int *nlon, int *nlat,
		       double *deltaLon, double *deltaLat);

int READTOPO_DEM_file(char *txt_inputfile, char *topoBathy_flg, 
		      double *rarray, int nlon, int nlat);


#endif
