/*
 * legendre.h
 */

void legendre_gauss_lobatto(int ngl, double *xgl, double *wgl);

void legendre_gauss(int ngl, double *xgl, double *wgl);

void legendre_basis(int ngl, double *xgl, double **psi, double **dpsi);

void lagrange_basis(int ngl, double *xgl, int nq, double *xnq,	\
		    double *wnq, double **psi, double **dpsi);

void legendre_poly(int n, double x, double *p0, double *p0_1, double *p0_2, double *p1, \
		   double *p1_1, double *p1_2, double *p2, double *p2_1, double *p2_2, \
		   double *p00, double *p00_1, double *p00_2);

