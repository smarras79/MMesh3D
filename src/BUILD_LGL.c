/***********************************************************************
 *This subroutine finds the Legendre-Gauss-Lobatto Roots.
 *
 * Original in Fortran 90 Written by Francis X. Giraldo on 7/08
 *           within the code NUMA3D 
 *           (Non hydrostatic Unified Model for the Atmosphere)
 *           Department of Applied Mathematics
 *           Naval Postgraduate School
 *           Monterey, CA 93943-5216
 *
 * C version by Simone Marras, August 2011
 *
 ***********************************************************************/
#include<stdio.h>
#include<math.h>

#include "BUILD_LGL.h"

void legendre_gauss_lobatto(int ngl, double *xgl, double *wgl)
{      
  
  int i,k, n, nh, kmax;
  double pi;
  double x, dx;
  double p00,p00_1,p00_2;
  double p0,p0_1,p0_2;
  double p1,p1_1,p1_2;
  double p2,p2_1,p2_2;

  pi=4.0*atan(1.0);
  n=ngl-1;
  nh=(n+1)/2;
  kmax = 20;
  
  //First Find Half of the Roots
  for (i=1; i<=nh; i++)
    {
      x=cos( (2.0*i - 1.0)/(2.0*n + 1.0)*pi );
      
      for( k=1; k<=kmax; k++)
	{
	  //Construct Legendre Polynomial and Derivatives
	  legendre_poly(n,&x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,	\
	  		&p2_2,&p00,&p00_1,&p00_2);
	  
	  //Get next Newton Iterative
	  dx=-(1.0-x*x)*p0_1/(-2.*x*p0_1 + (1.0-x*x)*p0_2);
	  
	  x=x + dx;
          if (abs(dx) < 1.0e-20)
	    continue;
	
	}//K loop

      xgl[n+2-i]=x;
      wgl[n+2-i]=2.0/( (double)(n*(n+1))*p0*p0 );
            
    }//I loop

  //Check for Zero
  if (n+1 != 2*nh)
    {
      x=0;
      legendre_poly(n,&x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,	\
		    &p2_2,&p00,&p00_1,&p00_2);
      
      xgl[nh+1]=x;
      wgl[nh+1]=2.0/( (double)(n*(n+1))*p0*p0 );
    }
  
  //Find Remainder of Roots via Symmetry
  for( i=1; i<=nh; i++)
    {
      xgl[i]= -xgl[n+2-i];
      wgl[i]=  wgl[n+2-i];
    }
  
  return;
}



/***********************************************************************
 *This subroutine finds the Legendre Polynomials and its Derivatives.
 *
 * Original in Fortran 90 Written by Francis X. Giraldo on 4/97
 *           within the code NUMA3D 
 *           (Non hydrostatic Unified Model for the Atmosphere)
 *           Department of Applied Mathematics
 *           Naval Postgraduate School
 *           Monterey, CA 93943-5216
 *
 * C version by Simone Marras, August 2011
 *
 ***********************************************************************/
void legendre_poly(int n, double *x,					\
		   double *p0, double *p0_1, double *p0_2,		\
		   double *p1, double *p1_1, double *p1_2,		\
		   double *p2, double *p2_1, double *p2_2,		\
		   double *p00, double *p00_1, double *p00_2)
{
 
  int j;
  double a,b,dx;
  
  //Construct Legendre Polynomials 00=N+1, 0=N, 1=N-1, 2=N-2, and derivs
  *p1=0;
  *p1_1=0;
  *p1_2=0;
    
  *p0=1;
  *p0_1=0;
  *p0_2=0;

  //Construct Nth Order Legendre Polynomial
  for( j=1; j<=n; j++)
    {
      
      *p2=*p1;
      *p2_1=*p1_1;
      *p2_2=*p1_2;

      *p1=*p0;
      *p1_1=*p0_1;
      *p1_2=*p0_2;

      a=(2.0*(double)(j)-1.0)/(double)(j);
      b=((double)(j)-1.0)/(double)(j);
      *p0=a* *x**p1 - b**p2;
      *p0_1=a*( *p1 + *x **p1_1 ) - b**p2_1;
      *p0_2=a*( 2.0**p1_1 + *x**p1_2 ) - b**p2_2;

      a=(2.0*(double)(j)+1.0)/((double)(j)+1.0);
      b=((double)(j))/((double)(j)+1.0);
      *p00=a* *x**p0 - b**p1;
      *p00_1=a*( *p0 + *x**p0_1 ) - b**p1_1;
      *p00_2=a*( 2.0**p0_1 + *x**p0_2 ) - b**p1_2;

      
    } //J loop
  
  return;
}
