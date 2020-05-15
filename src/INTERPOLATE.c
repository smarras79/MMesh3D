/***********************************************************
 *
 * INTERPOLATE.c 
 *
 * This file collects different functions used for the
 * interpolation of data points by means of an
 * interpolating surface.
 *
 * 1) LATTICE()
 *     Based on the algorithm BA, p. 231, of Lee et al.
 *     1997, IEEE Trans. Visualiz. Comput. Graphics, Vol 3
 *
 *
 ***********************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

void LATTICE(int m, int n, int nnodesx, int nnodesy, double **Q, float **PHI)
{
  /********************************************************
   * Input:
   * m,n: nodes in x and y for the lattice.
   * Q[1:npoints][1:3] = (xc,yc,zc): scattered data.
   * npoints: number of points of the scattered data.
   *
   * Output:
   * Lattice coordinates: PHI[1:m+3][1:n+3]
   ********************************************************/

  int i,j,k,l,ipoin;
  int npoints;

  float xc,yc,zc;
  float s,t;
  float delta[m+4][n+4],omega[m+4][n+4];
  float w[m+3][n+3];
  float B[4];               //Uniform cubic B-spline basis functions
  float phi[4][4];          //Local lattice of 16 poitns around point ipoin
  float sumWkl,sumWkl2;
  
  npoints = nnodesx*nnodesy;

  for(i=0; i<=m+3; i++){
    for(j=0; j<=n+3; j++){
      delta[i][j] = 0.0;
      omega[i][j] = 0.0;
    }
  }

  sumWkl = 0.0, sumWkl2 = 0.0;
  for(ipoin=1; ipoin<=npoints; ipoin++)
    {
      xc = Q[ipoin][2];
      yc = Q[ipoin][3];
      zc = Q[ipoin][4];
    
      i = floor(xc) - 1;
      j = floor(yc) - 1;
      
      s = xc - floor(xc);
      t = yc - floor(yc);
      k = (i + 1) - floor(xc);
      l = (j + 1) - floor(yc);
   
      B[0] = (1.0 - t)*(1.0 - t)*(1.0 - t)/6.0;
      B[1] = (3.0*t*t*t - 6.0*t*t + 4.0)/6.0;
      B[2] = (-3.0*t*t*t + 3.0*t*t + 3.0*t + 1.0)/6.0;
      B[3] = t*t*t/6.0;

      sumWkl2=0.0;
      for(k=0; k<=3; k++){
	for(l=0; l<=3; l++){ 
	  w[k][l] = B[k]*B[l];
	  sumWkl  = sumWkl + w[k][l];
	  sumWkl2 = sumWkl2 + sumWkl;
	}
      }
 
      for(k=0; k<=3; k++){
      	for(l=0; l<=3; l++){
	  phi[k][l] = w[k][l]*zc/(sumWkl2);
	  delta[k][l] = delta[k][l] + w[k][l]*w[k][l]*phi[k][l];
	  //omega[k+i][l+j] = omega[k+i][l+j] + w[k][l]*w[k][l];
	}
      }
   
         
    }//End ipoin
  
  for(i=0; i<=m+2; i++){
    for(j=0; j<=n+2; j++){
      if( omega[i][j] != 0.0)
	PHI[i][j] = delta[i][j]/(omega[i][j] + 1.0e-16);
      else
	PHI[i][j] = 0.0;
    }
  }
  
 for(i=0; i<=m+2; i++){
   for(j=0; j<=n+2; j++)
     printf(" %d %d %f", i,j,PHI[i][j]);
   printf("\n");
 }

  return;
}
