/*SURFACES.h*/

int SURFACES(int isurface, int nnodesu, int nnodesv, int nnodesw, int nelem, \
	     double **BDY_COORDS, char *problem[], double *parameters, double **COORDS, double *topo_surface,double deltaLon, double deltaLat);
/*
int SURFACES_HIGH_ORDER(int isurface, int nnodesu, int nnodesv, int nnodesw, int nop, int ngl, \
			int nelem, double **BDY_COORDS, char *problem[], double *parameters,   \
			double **COORDS, double *topo_surface,double deltaLon, double deltaLat);
*/
