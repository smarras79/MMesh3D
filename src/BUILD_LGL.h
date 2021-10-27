/***********************************************************************
 * BUILD_LGL.h
 ***********************************************************************/
int LegendreGaussNodesAndWeights(st_lgl lgl);
int LegendreGaussLobattoNodesAndWeights(st_lgl lgl);
st_legendre LegendreAndDerivative(int p, double x);
st_legendre LegendreAndDerivativeAndQ(int p, double x);
