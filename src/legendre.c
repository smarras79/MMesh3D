/*
 * legendre.c
 *
 * Contains:
 * legendre_gauss_lobatto()
 * legendre_gauss()
 * legendre_poly()
 * legendre_basis()
 * lagrange_basis()
 * lagrange_basis2()
 *
 */


#include"myinclude.h"
#include"global_vars.h"

/************************************************************************
 * 
 * This subroutine finds the Legendre-Gauss-Lobatto Roots.
 *
 * Fortran version: Written by F.X. Giraldo on 7/08 
 *                  Department of Applied Mathematics
 *                  Naval Postgraduate School
 *                  Monterey, CA 93943-5216
 *
 * C version: Written by S. Marras on 5/2013
 ************************************************************************/
void legendre_gauss_lobatto(int ngl, double *xgl, double *wgl)
{
  /*
   * Local variables:
   */
  int i,j,k;
  int kmax;
  int n, nh;
  
  double x;
  double dx;
  double p0;
  double p0_1;
  double p0_2;
  double p1;
  double p1_1;
  double p1_2;
  double p2;
  double p2_1;
  double p2_2;
  double p00;
  double p00_1;
  double p00_2;
  
  kmax = 20;
  n    = ngl - 1;
  nh   = 0.5*(n+1);
  dx   = 0.0; //initialize dx to be used before any computation
  
  for(i=1; i<=ngl; i++){
    xgl[i] = 0.0;
    wgl[i] = 0.0;
  }
  
  /*
   * First Find Half of the Roots
   */
  for(i=1; i<=nh; i++)
    {
      x = cos((2.0*i - 1.0)/(2.0*n + 1.0)*pi);
      
      for(k=1; k<=kmax; k++)
	{
	  /*
	   * Construct Legendre Polynomial and Derivatives
	   */
	  legendre_poly(n,x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,&p2_2,&p00,&p00_1,&p00_2);
	  
	  /*
	   * Get next Newton Iterative
	   */
	  dx = -(1.0 - x*x)*p0_1/(-2.0*x*p0_1 + (1.0 - x*x)*p0_2);
	  x  = x + dx;
	  if (fabs(dx) >= 1.0e-20 )
	    continue;
	  
	}
      xgl[n+2-i] = x;
      wgl[n+2-i] = 2.0/( (double)(n*(n+1))*p0*p0 );
    }
  
  /*
   * Check for Zero
   */
  if (n+1 != 2*nh)
    {
      x = 0.0;
      legendre_poly(n,x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,&p2_2,&p00,&p00_1,&p00_2);
      
      xgl[nh+1] = x;
      wgl[nh+1] = 2.0/( (double)(n*(n+1))*p0*p0 );
    }

  /*
   * Find Remainder of Roots via Symmetry
   */
  for(i=1; i<=nh; i++)
    {
      xgl[i] = -xgl[n+2-i];
      wgl[i] =  wgl[n+2-i];
    }
  
  return;
}

/************************************************************************
 *
 * This subroutine finds the Legendre-Gauss Roots.
 * Fortran version: Written by F.X. Giraldo on 9/97 
 *           Naval Research Laboratory
 *           Monterey, CA 93943
 *
 * C version: Written by S. Marras on 5/2013
 *
 ************************************************************************/
void legendre_gauss(int ngl, double *xgl, double *wgl)
{
 
  /*
   * Local variables:
   */
  int i,j,k;
  int kmax;
  int n, nh;
  
  double x,dx;
  double p0;
  double p0_1;
  double p0_2;
  double p1;
  double p1_1;
  double p1_2;
  double p2;
  double p2_1;
  double p2_2;
  double p00;
  double p00_1;
  double p00_2;
  
  kmax = 20;
  n    = ngl-1;
  nh   = 0.5*(n+1);
  dx   = 0.0; //initialize dx to be used before any computation

  /*
   * First Find Half of the Roots
   */
  for(i=1; i<=nh; i++)
    {
      x = cos( (2.*i - 1.)/(2.*n + 1.)*pi );

      for(k=1; k<=kmax; k++)
        {
	  
	  /*
	   * Construct Legendre Polynomial and Derivatives
	   */
	  legendre_poly(n,x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,&p2_2,&p00,&p00_1,&p00_2);
	  
	  /*
	   * Get next Newton Iterative
	   */
	  dx = -p00/p00_1;
	  x  = x + dx;
	  if (fabs(dx) >= 1.0e-20)
	    continue;
	}
      
      xgl[n+2-i] = x;
      wgl[n+2-i] = 2.0/( (1.0-x*x)*p00_1*p00_1 );
    }

  /*
   * Check for Zero
   */
  if (n+1 != 2*nh)
    {
      x = 0.0;
      legendre_poly(n,x,&p0,&p0_1,&p0_2,&p1,&p1_1,&p1_2,&p2,&p2_1,&p2_2,&p00,&p00_1,&p00_2);

      xgl[nh+1] = x;
      wgl[nh+1] = 2.0/( (1.0 - x*x)*p00_1*p00_1 );
    }

  /*
   * Find Remainder of Roots via Symmetry
   */
  for(i=1; i<=nh; i++)
    {
      xgl[i] = -xgl[n+2-i];
      wgl[i] =  wgl[n+2-i];
    }
   
  return;
}


/************************************************************************
 * 
 * This subroutine finds the Legendre Cardinal Basis Functions and 
 * their Derivatives.
 * Fortran version: Written by F.X. Giraldo on 4/97 
 *           Naval Research Laboratory
 *           Monterey, CA 93943
 *
 * C version: Written by S. Marras on 5/2013
 ************************************************************************/
void legendre_basis(int ngl, double *xgl, double **psi, double **dpsi)
{
  
  /*
   * Local variables:
   */
  int i,j,n;
  double ksi;
  double xj;

  double p0i;
  double p0i_1;
  double p0i_2;
  double p1i;
  double p1i_1;
  double p1i_2;
  double p2i;
  double p2i_1;
  double p2i_2;

  double p0j;
  double p0j_1;
  double p0j_2;
  double p1j;
  double p1j_1;
  double p1j_2;
  double p2j;
  double p2j_1;
  double p2j_2;
  double p00;
  double p00_1;
  double p00_2;
  
  n = ngl-1;

  for(j=1; j<=ngl; j++)
    {
      xj = xgl[j];
      legendre_poly(n,xj,&p0j,&p0j_1,&p0j_2,&p1j,&p1j_1,&p1j_2,&p2j,&p2j_1,&p2j_2,&p00,&p00_1,&p00_2);
        
      for(i=1; i<=ngl; i++)
	{
	  ksi = xgl[i];
	  legendre_poly(n,ksi,&p0i,&p0i_1,&p0i_2,&p1i,&p1i_1,&p1i_2,&p2i,&p2i_1,&p2i_2,&p00,&p00_1,&p00_2);
	  if (i == j)
	    psi[i][j] = 1.0;
	  else
	    psi[i][j] = 0.0;
            
	  if (i == j)
	    {
	      if (i != 1 && i != ngl)
		dpsi[i][j] = 0;		
	      else if (i == 1) 
		dpsi[i][j] = -(double)(n*(n+1))/4.0;
	      else if (i == ngl)
		dpsi[i][j] = (double)(n*(n+1))/4.0;
	    }
	  else if (i != j)
	    dpsi[i][j] = p0j/(p0i*(xj-ksi));
	  
	}
    }

  return;
}


/************************************************************************
 * 
 * This subroutine finds the Legendre Cardinal Basis Functions and 
 * their Derivatives in Lagrange Polynomial Form.
 * Fortran version: Written by F.X. Giraldo on 4/97 
 *           Naval Research Laboratory
 *           Monterey, CA 93943
 *
 * C version: Written by S. Marras on 5/2013
 ************************************************************************/
void lagrange_basis(int ngl, double *xgl, int nq, double *xnq,	\
		    double *wnq, double **psi, double **dpsi)
{
  /*
   * Local variables:
   */
  int i,j,k,l;
  
  double xi,xj,xk,xl;
  double ksi;
  double ddpsi;

  /*
   * Get High Order Roots
   */
  legendre_gauss_lobatto(nq, xnq, wnq);

  /*
   * Do Quadrature
   */    
  for(l=1; l<=nq; l++)
    {
      xl = xnq[l];

      /*
       * Construct Bases
       */
      for(i=1; i<=ngl; i++)
	{
	  ksi       = xgl[i];
	  psi[i][l] = 1.0;
	  dpsi[i][l]= 0.0;

	  for(j=1; j<=ngl; j++)
	    {
	      xj = xgl[j];

	      /*
	       * Basis Functions
	       */
	      if (j != i)
		psi[i][l] = psi[i][l]*( xl - xj )/( ksi - xj );
		
	      ddpsi = 1.0;
               
	      /*
	       * Derivative of Basis Functions
	       */
	      if (j != i)
		{
		  for(k=1; k<=ngl; k++)
		    {
		      xk = xgl[k];
			
		      if (k != i && k != j)
			ddpsi = ddpsi*( xl - xk )/( ksi - xk );
			
		    }
		  dpsi[i][l] = dpsi[i][l] + ddpsi/( ksi - xj );
		}//end if
   
	    }//end for j
	}//end for i
    }//end for l
   
  return;
}


/************************************************************************
 *
 * This subroutine finds the Legendre Polynomials and its Derivatives.
 * Fortran version: Written by F.X. Giraldo on 4/97 
 *          Naval Research Laboratory
 *          Monterey, CA 93943
 *
 * C version: written by S. Marras on 5/2013
 *
 ************************************************************************/
void legendre_poly(int n, double x, double *p0, double *p0_1, double *p0_2, double *p1, \
		   double *p1_1, double *p1_2, double *p2, double *p2_1, double *p2_2, \
		   double *p00, double *p00_1, double *p00_2)
{

  /*
   * Local variables:
   */
  int j;
  double a,b;

  /*
   * Construct Legendre Polynomials 00=N+1, 0=N, 1=N-1, 2=N-2, and derivs
   */
  *p1   = 0.0;
  *p1_1 = 0.0;
  *p1_2 = 0.0;

  *p0   = 1.0;
  *p0_1 = 0.0;
  *p0_2 = 0.0;

  /*
   * Construct Nth Order Legendre Polynomial
   */
  for(j=1; j<=n; j++)
    {
      *p2    = *p1;
      *p2_1  = *p1_1;
      *p2_2  = *p1_2;

      *p1    = *p0;
      *p1_1  = *p0_1;
      *p1_2  = *p0_2;

      a     = (2.0*(double)(j) - 1.0)/(double)(j);
      b     = ((double)(j) - 1.0)/(double)(j);
      
      *p0    = a*x**p1 - b**p2;
      *p0_1  = a*( *p1 + x**p1_1 ) - b**p2_1;
      *p0_2  = a*( 2.0**p1_1 + x**p1_2 ) - b**p2_2;

      a     = (2.0*(double)(j) + 1.0)/((double)(j) + 1.0);
      b     = ((double)(j))/((double)(j) + 1.0);
      *p00   = a*x**p0 - b**p1;
      *p00_1 = a*( *p0 + x**p0_1 ) - b**p1_1;
      *p00_2 = a*( 2.0**p0_1 + x**p0_2 ) - b**p1_2;

    }

  return;
}
