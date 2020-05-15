int BUILD_GRID_SPHERE(int rank);
int init_hex(double **xcoord, double **ycoord, int **intma, int nelem, int nel, int nx, int ny, int ngl);
int hex_quad(double **coord, int **intma, int nel, int ngl, int nface, int nelem, \
	     int npoin, int nelem0, int npoin0, int nx, int ny);
int create_grid_sphere_hex(double **coord, int **intma, int *intma_r, \
			   int ***intma_s, int nel, int nelem_r, int nelem_s, \
			   int npoin_r, int npoin_s, double rmin, double rmax);
int create1d_grid(double *r, int *intmar, int nelr, int npoinr, double rmin, double rmax);
