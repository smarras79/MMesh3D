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
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>

#include "mystructs.h"
#include "global_vars.h"
#include "BUILD_LGL.h"

#define TOL 4*DBL_EPSILON
#define NITER 20


/***********************************************************************
 * Evaluate by recursion, the Legendre polynomial of order p
 * and its Derivatives at coordinate x
 *
 * L_{p}  --> legendre of order p
 * L'_{p} --> dlegendr of order p
 * 
 * Algorithm 22 of Kopriva's book
 *
 * Simone Marras, October 2021
 *
 ***********************************************************************/   
st_legendre LegendreAndDerivative(int p, double x)
{
    st_legendre Legendre;
    
    int i, k;
    double a, b;
    double Lm1, dLm1, L0, dL0, L, dL;

    Lm1  = 0;
    dLm1 = 0;
    
    L0  = 1;
    dL0 = 0;
    
    //Construct Nth Order Legendre Polynomial
    for (k = 1; k<=p; k++){
      
	a = (double)(2*k - 1)/k;
	b = (double)(k - 1)/k;
      
	L  = a*x*L0 - b*Lm1;
	dL = a*( L0 + x*dL0 ) - b*dLm1;

	Lm1 = L0;
	L0  = L;

	dLm1 = dL0;
	dL0  = dL;
    }

    Legendre.legendre  = L;
    Legendre.dlegendre = dL;
  
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
st_legendre LegendreAndDerivativeAndQ(int p, double x)
{
    st_legendre Legendre;
    
    int    i, k;
    double a, b;
    double  Lm2, dLm2, Lm1, dLm1, Lp1, dLp1, L, dL;

    k = 2;
    
    Lm2  = 1.0;
    Lm1  = x;
    
    dLm2 = 0.0;
    dLm1 = 1.0;
        
    //Construct Nth Order Legendre Polynomial
    for (k = 2; k<=p; k++){
      
	a = (double)(2*k - 1)/k;
	b = (double)(k - 1)/k;
      
	L  = a*x*Lm1 - b*Lm2;
	dL = dLm2 + (2*k - 1)*Lm1;

	Lm2 = Lm1;
	Lm1 = L;

	dLm2 = dLm1;
	dLm1 = dL;
    }
    k = p + 1;

    Lp1  = a*x*L - b*Lm1;
    dLp1 = dLm1 + (2*k - 1)*L;

    Legendre.legendre  = L;
    Legendre.dlegendre = dL;
    
    Legendre.q  = Lp1  - Lm1;
    Legendre.dq = dLp1 - dLm1;
    
    return Legendre;
}

/***********************************************************************
 * Calculate the Legendre-Gauss nodes and weights
 *
 * Algorithm 23 of Kopriva's book
 *
 * Simone Marras, October 2021
 *
 ***********************************************************************/
int LegendreGaussNodesAndWeights(st_lgl lgl)
{
    st_legendre Legendre;
    
    int p;
    double x0, x1;
    double w0, w1;
    double xj;
    double Delta;
    
    //Polynomial order
    p = nop;
    
    printf(" NOP = %d\n", p);
    
    if (p == 0) {
	x0 = 0.0;
	w0 = 2.0;
    } else if (p == 1) {
	x0 = -sqrt(1.0/3.0);
	w0 = 1.0;
	x1 = -x0;
	w1 =  w0;
    } else {
	for (int j=0; j<=((p + 1)/2) - 1; j++)
	    {
		
		xj = -cos((2*j + 1)/(2*p + 2)*PI);
		lgl.coords[j] = xj;
		
		for (int k=0; k<=NITER; k++)
		    {
			Legendre =LegendreAndDerivative(p + 1, xj);
			Delta         = -Legendre.legendre/Legendre.dlegendre;
			xj            = xj + Delta;
			if (fabs(Delta) <= TOL*fabs(xj)) break;
		    }
		lgl.coords[j] = xj;
		Legendre =LegendreAndDerivative(p + 1, xj);
		lgl.coords[p - j] = -xj;
		
		double xj2 = xj*xj;
		double dL2 = Legendre.dlegendre*Legendre.dlegendre;
		lgl.weights[j]     = 2/((1 - xj2)*dL2);
		lgl.weights[p - j] = lgl.weights[j];
		
	    }
    }

    if ((p % 2) == 0){
	Legendre = LegendreAndDerivative(p + 1, 0.0);
	lgl.coords[p/2] = 0.0;

	double dL2 = Legendre.dlegendre*Legendre.dlegendre;
	lgl.weights[p/2] = 2.0/dL2;
    }

    for (int j=0; j<p+1; j++){
	printf("TOL = %.16e. X,W: %d, %.8f %.8f\n", TOL, j, lgl.coords[j], lgl.weights[j]);
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
 *
 ***********************************************************************/
int LegendreGaussLobattoNodesAndWeights(st_lgl lgl)
{
    st_legendre Legendre;
    
    int p;
    double x0, x1, xP;
    double w0, w1, wP;
    double xj;
    double Delta;

    p = nop;

    if (p == 1) {
	x0 = -1.0;
	w0 =  1.0;
	x1 =  1.0;
	w1 =  w0;
    } else {	
	x0 = -1.0;
	w0 =  2.0/(p*(p+1));
	xP =  1.0;
	wP =  w0;
	
	for (int j=1; j<=((p + 1)/2) - 1; j++)
	    {		
		xj = -cos((j + 0.25)*PI/p - 3.0/(8*p*PI*(j + 0.25)));
		lgl.coords[j] = xj;
		
		for (int k=0; k<=NITER; k++)
		    {
			Legendre = LegendreAndDerivativeAndQ(p, xj);
			Delta         = Legendre.q/Legendre.dq;
			xj            = xj + Delta;
			if (fabs(Delta) <= TOL*fabs(xj)) break;
		    }
		lgl.coords[j] = xj;
		Legendre = LegendreAndDerivativeAndQ(p, xj);
		lgl.coords[p - j] = -xj;
		
		double xj2 = xj*xj;
		double L2  = Legendre.legendre*Legendre.legendre;
		lgl.weights[j]     = 2/(p*(p + 1)*L2);
		lgl.weights[p - j] = lgl.weights[j];
	    }
    }
	
    if ((p % 2) == 0){
	Legendre        = LegendreAndDerivativeAndQ(p, 0.0);
	lgl.coords[p/2] = 0.0;
	    
	double L2        = Legendre.legendre*Legendre.legendre;
	lgl.weights[p/2] = 2/(p*(p + 1)*L2);
    }
	
    for (int j=0; j<lgl.size+1; j++){
	printf("TOL = %.16e. X,W: %d, %.8f %.8f\n", TOL, j, lgl.coords[j], lgl.weights[j]);
    }
}

