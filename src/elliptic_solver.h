void elliptic_solver(double **COORDS, int **CONN, int *ELTYPE,		\
		     double **x, double **y, double **z, double **BDYFLAG, \
		     double **bottomSide, double **rightSide, double **topSide, double **leftSide, \
		     double *bottomSide_der, double *rightSide_der, double *topSide_der, double *leftSide_der, \
		     int nnodesx, int nnodesy, int nnodesz);

void gauss_seidel_xy( double **x, double **y, double **z, int nnodesx, int nnodesy, int nnodesz, double *err);

void compute_f( double **x, double **y, int nnodesx, int nnodesy, double **f);

void gauss_seidel_TTM(int iter, double ***x, double ***y, double ***z, double err);

