#ifndef LINSPACE_H
#define LINSPACE_H

void flinspace(float vmin, float vmax, int n, float *VECT, char *force_flag);

void ilinspace(int vmin, int vmax, int n, int *VECT);

void dlinspace(double vmin, double vmax, int n, double *VECT, char *force_flag);

void linspace_test(float vmin, float vmax, int n, float *VECT);

#endif
