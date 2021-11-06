#ifndef BUILD_LGL_H
#define BUILD_LGL_H

/***********************************************************************
 * BUILD_LGL.h
 ***********************************************************************/
st_lgl BUILD_LGL(size_t p);

st_legendre LegendreAndDerivative(size_t p, double x);
st_legendre LegendreAndDerivativeAndQ(size_t p, double x);

int LegendreGaussNodesAndWeights(st_lgl lgl, size_t p);
int LegendreGaussLobattoNodesAndWeights(st_lgl lgl, size_t p);
int BarycentricWeights(st_lgl lgl, size_t p);

#endif
