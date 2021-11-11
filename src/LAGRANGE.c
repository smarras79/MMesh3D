/***********************************************************************
 *
 * functions to compute the barycentric weights and Lagrange interpolation
 * 
 * Simone Marras, October 2021
 *
 * References:
 * [1] D. Kopriva "Implementing spectral methods for PDEs" Springer
 *
 ***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

//User defined structures
#include "MYSTRUCTS.h"

//Function headers
#include "LAGRANGE.h"
#include "PRINT.h"

//Local constants
#define TOL 4*DBL_EPSILON
#define NITER 100

/***********************************************************************
 * BarycentricWeights
 * 
 * Algorithm 30 of Kopriva's book
 *
 * Simone Marras, October 2021
 ***********************************************************************/
int BarycentriWeights(size_t p, st_lgl lgl)
{
    for (int j = 0; j < p;  j++)
	lgl.weights[j] = 1.0;


    VIEW_dVECT("LGL.WEIGHTS", lgl.weights, 0, p-1);
    
    return 0; 
}
