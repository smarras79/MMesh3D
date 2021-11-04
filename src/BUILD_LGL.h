/***********************************************************************
 * BUILD_LGL.h
 ***********************************************************************/
int BUILD_LGL(size_t p, st_lgl lgl);

int LegendreGaussNodesAndWeights(st_lgl lgl, size_t p);

int LegendreGaussLobattoNodesAndWeights(st_lgl lgl, size_t p);

st_legendre LegendreAndDerivative(size_t p, double x);
