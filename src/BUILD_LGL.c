/***********************************************************************
 *
 * functions to compute the (L)GL nodes and weights for spectral elements
 * 
 * Simone Marras, August 2011
 *
 * References:
 * [1] D. Kopriva "Implementing spectral methods for PDEs" Springer
 * [2] F. X. Giraldo "An introduction to element0-based Galerkin ..." Springer
 *
 ***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

//User defined structures
#include "MYSTRUCTS.h"
#include "MESH.h"

//Global variable declarations
#include "GLOBAL_VARS.h"

//Function headers
#include "ALMOST_EQUAL.h"
#include "BUILD_LGL.h"
#include "PRINT.h"

//Local constants
#define TOL 4*DBL_EPSILON
#define NITER 100


/***********************************************************************
 * Evaluate by recursion, the Legendre polynomial of order p
 * and its Derivatives at coordinate x
 *
 * L_{p}  --> legendre of order p
 * L'_{p} --> dlegendr of order p
 * 
 * and
 *
 * q  = L_{p+1}  -  L_{p-1}
 * q' = L'_{p+1} -  L'_{p-1} 
 * 
 * Algorithm 22+24 of Kopriva's book
 *
 * Simone Marras, October 2021
 *OK
 ***********************************************************************/
st_lgl BUILD_LGL(size_t p)
{
    int ivalue;
    st_legendre Legendre;
    st_lgl lgl;
		
    lgl.size    = ngl;
    lgl.ksi     = (double*) malloc( sizeof(double) * lgl.size);
    lgl.weights = (double*) malloc( sizeof(double) * lgl.size);
    
    //LGL nodes
    ivalue = LegendreGaussLobattoNodesAndWeights(lgl, nop);
    
    //LG nodes
    //ivalue = LegendreGaussNodesAndWeights(lgl, nop);
    
    return lgl;
}

/***********************************************************************
 * Evaluate by recursion, the Legendre polynomial of order p
 * and its Derivatives at coordinate x
 *
 * L_{p}  --> legendre of order p
 * L'_{p} --> dlegendr of order p
 * 
 * and
 *
 * q  = L_{p+1}  -  L_{p-1}
 * q' = L'_{p+1} -  L'_{p-1} 
 * 
 * Algorithm 22+24 of Kopriva's book
 *
 * Simone Marras, October 2021
 *OK
 ***********************************************************************/
st_legendre LegendreAndDerivative(size_t p, double x)
{
    st_legendre Legendre;
    
    int i, k;
    double a, b;
    double Lm2, dLm2, Lm1, dLm1, L0, dL0, L, dL, Lp1, dLp1;
    
    if (p == 0)
	{
	    L  = 1.0;
	    dL = 0.0;
	}
    else if (p == 1)
	{
	    L  = x;
	    dL = 1.0;
	}
    else
	{
	    Lm2  = 1.0;
	    Lm1  = x;
	    dLm2 = 0.0;
	    dLm1 = 1.0;
	    
	    //Construct Nth Order Legendre Polynomial
	    for (k = 2; k<=p; k++){
      
		a = (double)(2.0*k - 1.0)/k;
		b = (double)(k - 1.0)/k;
      
		L  = a*x*Lm1 - b*Lm2;
		dL = dLm2 + (2.0*k - 1.0)*Lm1;

		Lp1  = a*x*L - b*Lm1;
		dLp1 = dLm1 + (2.0*k - 1.0)*L;
    
		Lm2 = Lm1;
		Lm1 = L;

		dLm2 = dLm1;
		dLm1 = dL;
	    }
	}
    
    k = p + 1;
    
    a = (double)(2.0*k - 1.0)/k;
    b = (double)(k - 1.0)/k;
    
    Legendre.legendre  = L;
    Legendre.dlegendre = dL;
    
    Legendre.q  =  Lp1 -  Lm2;
    Legendre.dq = dLp1 - dLm2;
  
    return Legendre;
}

/***********************************************************************
 * Evaluate by recursion, the Legendre polynomial of order p
 * its Derivatives at coordinate x, and q:
 *
 * L_{p}  --> legendre of order p
 * L'_{p} --> dlegendr of order p
 * q  = L_{p+1}  -  L_{p-1}
 * q' = L'_{p+1} -  L'_{p-1} 
 *
 * Algorithm 24 of Kopriva's book
 *
 * Simone Marras, October 2021
 *
 ***********************************************************************/
st_legendre LegendreAndDerivativeAndQ(size_t p, double x)
{
    st_legendre Legendre;
    
    int i, k;
    double a, b;
    double Lm2, dLm2, Lm1, dLm1, L0, dL0, L, dL, Lp1, dLp1;
    
    if (p == 0)
	{
	    L  = 1.0;
	    dL = 0.0;
	}
    else if (p == 1)
	{
	    L  = x;
	    dL = 1.0;
	}
    else
	{
	    Lm2  = 1.0;
	    Lm1  = x;
	    dLm2 = 0.0;
	    dLm1 = 1.0;
	    
	    //Construct Nth Order Legendre Polynomial
	    for (k = 2; k<=p; k++){
      
		a = (double)(2.0*k - 1.0)/k;
		b = (double)(k - 1.0)/k;
      
		L  = a*x*Lm1 - b*Lm2;
		dL = dLm2 + (2.0*k - 1.0)*Lm1;
		
		a = (double)(2.0*(k+1) - 1.0)/(k+1);
		b = (double)(k+1 - 1.0)/(k+1);
		Lp1  = a*x*L - b*Lm1;
		dLp1 = dLm1 + (2.0*k - 1.0)*L;
		
		Lm2 = Lm1;
		Lm1 = L;

		dLm2 = dLm1;
		dLm1 = dL;
	    }
	}
    
    Legendre.legendre  = L;
    Legendre.dlegendre = dL;
    
    Legendre.q  =  Lp1 -  Lm2;
    Legendre.dq = dLp1 - dLm2;
		
    return Legendre;
}

/***********************************************************************
 * Calculate the Legendre-Gauss nodes and weights
 *
 * Algorithm 23 of Kopriva's book
 *
 * Simone Marras, October 2021
 * OK
 ***********************************************************************/
int LegendreGaussNodesAndWeights(st_lgl lgl, size_t p)
{
    st_legendre Legendre;
    
    double x0, x1;
    double w0, w1;
    double xj, wj;
    double Delta;
    double xj2, dL2;
    
    if (p == 0) {
	x0 = 0.0;
	w0 = 2.0;
	lgl.ksi[p]  = x0;
	lgl.weights[p] = w0;
    } else if (p == 1) {
	x0 = -sqrt(1.0/3.0);
	w0 = 1.0;
	x1 = -x0;
	w1 =  w0;
	lgl.ksi[0]  = x0;
	lgl.weights[0] = w0;
	lgl.ksi[p]  = x1;
	lgl.weights[p] = w1;
    } else {
	
	for (int j=0; j<(int)(p + 1)/2 + 1; j++)
	    {
		xj = -cos(PI*(2.0*j + 1.0)/(2.0*p + 2.0));
		lgl.ksi[j] = xj;
		
		for (int k=0; k<=NITER; k++)
		    {
			Legendre =  LegendreAndDerivative(p+1, xj);
			Delta    = -Legendre.legendre/Legendre.dlegendre;
			xj       = xj + Delta;
		
			if (fabs(Delta) <= TOL*fabs(xj)) break;
		    }
		Legendre      = LegendreAndDerivative(p + 1, xj);
			
		lgl.ksi[j]      =  xj;
		lgl.ksi[p - j]  = -xj;
		xj2                =  xj*xj;
		dL2                = Legendre.dlegendre*Legendre.dlegendre;
		lgl.weights[j]     = 2.0/((1 - xj2)*dL2);
		lgl.weights[p - j] = lgl.weights[j];
	    }
    }
    
    if ((p % 2) == 0){
	Legendre = LegendreAndDerivative(p + 1, 0.0);
	lgl.ksi[p/2] = 0.0;

	dL2 = Legendre.dlegendre*Legendre.dlegendre;
	lgl.weights[p/2] = 2.0/dL2;
    }
    
    
    for (int j=0; j<=p; j++){
	printf(" # LG nodes: X_LG[%d] = %.16f, w_LG[%d] = %.16f\n", j, lgl.ksi[j], j, lgl.weights[j]);
    }
    
    return 0;
}


/***********************************************************************
 * Calculate the Legendre-Gauss-Lobatto nodes and weights as the roots 
 * of Legendre polynomial derivatives of order p
 * 
 * Algorithm 25 of Kopriva's book
 *
 * Simone Marras, October 2021
 *OK
 ***********************************************************************/
int LegendreGaussLobattoNodesAndWeights(st_lgl lgl, size_t p)
{
    st_legendre Legendre;
    
    double x0, x1, xP;
    double w0, w1, wP;
    double xj, xj2, L2;
    double Delta;

    printf(" # Compute LGL nodes ... \n");
    
    for (int j=0; j<=p; j++) {
	lgl.ksi[j]  = 0.0;
	lgl.weights[j] = 2.0;
    }
    
    if (p == 1) {
	x0 = -1.0;
	w0 =  1.0;
	x1 =  1.0;
	w1 =  w0;
	lgl.ksi[0]  = x0;
	lgl.ksi[p]  = x1;
	lgl.weights[0] = w0;
	lgl.weights[p] = w1;
    } else {	
	x0 = -1.0;
	w0 =  (double)2.0/(p*(p + 1));
	xP =  1.0;
	wP =  w0;
	lgl.ksi[0]  = x0;
	lgl.ksi[p]  = xP;
	lgl.weights[0] = w0;
	lgl.weights[p] = wP;
	
	for (int j=1; j<(int)(p + 1)/2 + 1; j++)
	    {		
		xj            = -cos((j + 0.25)*PI/p - 3.0/(8.0*p*PI*(j + 0.25)));
		lgl.ksi[j] = xj;
		
		for (int k=0; k<=NITER; k++)
		    {
			Legendre      =  LegendreAndDerivativeAndQ(p, xj);
			Delta         = -Legendre.q/Legendre.dq;
			xj            = xj + Delta;
		
			if (fabs(Delta) <= TOL) break;
		    }
		Legendre           = LegendreAndDerivativeAndQ(p, xj);
		lgl.ksi[j]      =  xj;
		lgl.ksi[p - j]  = -xj;
		xj2                =  xj*xj;
		L2                 = Legendre.legendre*Legendre.legendre;
		lgl.weights[j]     = 2/(p*(p + 1)*L2);
		lgl.weights[p - j] = lgl.weights[j];
	    }
    }
    
    if ((p % 2) == 0){
	Legendre        = LegendreAndDerivative(p, 0.0);
	lgl.ksi[p/2] = 0.0;
	    
	L2        = Legendre.legendre*Legendre.legendre;
	lgl.weights[p/2] = 2/(p*(p + 1)*L2);
    }

    for (int j=0; j<=p; j++){
	printf(" #\t LGL nodes: X_LG[%d] = %.16f; w_LG[%d] = %.16f\n", j, lgl.ksi[j], j, lgl.weights[j]);
    }

    printf(" # Compute LGL nodes           ... DONE\n");
    return 0;
}


/***********************************************************************
 * BarycentricWeights
 * 
 * Algorithm 30 of Kopriva's book
 *
 * Simone Marras, October 2021
 ***********************************************************************/
int BarycentricWeights(st_lgl lgl, size_t p)
{
    //VIEW_dVECT("LGL.WEIGHTS", lgl.weights, 0, lgl.size-1);
    
    printf(" # Compute barycentric weights ... \n");
    for (int j = 0; j <= p;  j++)
	lgl.weights[j] = 1.0;

    for (int j = 1; j <= p;  j++) {
	for (int k = 0; k <= j - 1;  k++) {
	    double xk = lgl.ksi[k];
	    double xj = lgl.ksi[j];

	    double wk = lgl.weights[k];
	    double wj = lgl.weights[j];
	    
	    wk = wk*(xk - xj);
	    wj = wj*(xj - xk);
	    
	    lgl.weights[k] = wk;
	    lgl.weights[j] = wj;
	}	
    }
    
    for (int j = 0; j <= p;  j++)
	lgl.weights[j] = 1.0/lgl.weights[j];
    
    for (int j=0; j<=p; j++){
	printf(" #\t Barycentryc weights: w_LGL[%d] = %.16f\n", j, lgl.weights[j]);
    }
    printf(" # Compute barycentric weights ... DONE\n");
    return 0; 
}
