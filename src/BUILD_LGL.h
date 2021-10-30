/***********************************************************************
 * BUILD_LGL.h
 ***********************************************************************/
int LegendreGaussNodesAndWeights(st_lgl lgl, size_t p);
int LegendreGaussLobattoNodesAndWeights(st_lgl lgl, size_t p);
st_legendre LegendreAndDerivative(int p, double x);
st_legendre LegendreAndDerivativeAndQ(int p, double x);
int LegendreGaussNodesAndWeights_giraldo(st_lgl lgl, size_t p);
int LegendreGaussNodesAndWeights_ks_appendixB2(st_lgl lgl, size_t p);
