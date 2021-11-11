#ifndef INTERPOLATE_H
#define INTERPOLATE_H

double LagrangeInterpolation(double x, st_lgl lgl, double *f, size_t p);

void LATTICE(int m, int n, int nnodesx, int nnodesy, double **Q, float **PHI);

#endif
