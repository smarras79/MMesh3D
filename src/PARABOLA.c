// parabola.c

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//User-defined structs and variables
#include "GAUSSJ.h"
#include "GLOBAL_VARS.h" //GLOBAL VARIABLES
#include "MYDEFINE.h"    //GLOBAL CONSTANTS
#include "NRUTIL.h"
#include "PARABOLA.h"

void parabola(double *point1, double *point2, double *point3, char *function_type, double *a, double *b, double *c)
{
  /* Function that computes the coefficients a,b, and c of a parabola
   * given by the generic formula:
   *
   *	z(ordinate) = a*ordinate^2 + b*ordinate + c
   *	(e.g. y(x) = a*x^2 + b*x * c where x could be y or any other coordinate depending on your specific system of ref.)
   *
   * where "ordinate" is known given three passing points. 
   *
   *
   *	INPUTS:
   *	point1, point2, point3 are arrays of size [2]X[1] with containing the X,and Y of the point.
   *	function_type = "y(x)" indicating that we are using the canonical y = a*x^2 + b*x + c
   *	function_type = "x(y)" indicating that we are using the flipped   x = a*y^2 + b*y + c (e.g. parabola along the x axis instead of vertical)
   *
   * Given three known points (point)_{1,2,3} defined as point_1 = (absissa1, ordinate1), 
   * if parabola z(absissa) = a*absissa^2 + b*absissa + c
   * passes thorugh them, the coefficients a,b,c are obtained by solving the system:
   * 
   *	[a b c]' = inv(A)*[ordinate1 ordinate2 ordinate3]'
   * where A is:
   *	    [absissa1^2 absissa1 1]
   *	A = [absissa2^2 absissa2 1]
   *	    [absissa3^2 absissa3 1]
   *
   *
   * Simone Marras
   */
  
  int i;
  
  double tmp;
  double **A, uknown[3], **Z;
  
  double abscissa1, abscissa2, abscissa3, ordinate1, ordinate2, ordinate3;
  
  A = dmatrix(1,3,1,3);
  Z = dmatrix(1,3,1,1);
  
  //Initialize values:
  A[1][1]=0.0, A[1][2]=0.0, A[1][3]=0.0;
  A[2][1]=0.0, A[2][2]=0.0, A[2][3]=0.0;
  A[3][1]=0.0, A[3][2]=0.0, A[3][3]=0.0;
  
  Z[1][1]=0.0;
  Z[2][1]=0.0;
  Z[3][1]=0.0;
  
  uknown[0]=0.0;
  uknown[1]=0.0;
  uknown[2]=0.0;
  //End INITIALIZATION.
  
  //Switching of the coordinates order because the parabola is x(y) = ay^2 + by + c
  if( !strcmp(function_type, "x(y)"))
    {
      tmp = point1[1];
      point1[1] = point1[2];
      point1[2] = tmp;
      tmp = point2[1];
      point2[1] = point2[2];
      point2[2] = tmp;
      tmp = point3[1];
      point3[1] = point3[2];
      point3[2] = tmp;
    }
  
  abscissa1 = point1[1];
  ordinate1 = point1[2];
  
  abscissa2 = point2[1];
  ordinate2 = point2[2];
  
  abscissa3 = point3[1];
  ordinate3 = point3[2];
  
  A[1][1] = abscissa1*abscissa1;
  A[1][2] = abscissa1;
  A[1][3] = 1.0;
  
  A[2][1] = abscissa2*abscissa2;
  A[2][2] = abscissa2;
  A[2][3] = 1.0;
  
  A[3][1] = abscissa3*abscissa3;
  A[3][2] = abscissa3;
  A[3][3] = 1.0;
  
  Z[1][1] = ordinate1;
  Z[2][1] = ordinate2;
  Z[3][1] = ordinate3;
  
  dgaussj(A,3,Z,1);
  
  *a = Z[1][1];
  *b = Z[2][1];
  *c = Z[3][1];
  
  free_dmatrix(A,1,3,1,3);
  free_dmatrix(Z,1,3,1,1);
  
  return;
}


/***********************************************************************************
 *
 * parabola2pts.c computes a,b,c as parabola.c, but when the mid point is not
 *                passed as point2 because it was not defined. In this case, 
 *                the mid point is computed inside this function as the medium value
 *                of point1 and point2
 *
 **********************************************************************************/
void parabola2pts(double *point1, double *point3, char *function_type, \
		  double *a, double *b, double *c)
{

	
	int i;
	
	double tmp;
	double **A, uknown[3], **Z;
	
	double *point2;

	double abscissa1, abscissa2, abscissa3, ordinate1, ordinate2, ordinate3;
	
	A = dmatrix(1,3,1,3);
	Z = dmatrix(1,3,1,1);
	point2 = dvector(1,3); //(x,y,z)
	
	//Compute the inexistent point 2:
	point2[1] = (point1[1] + point3[1])/2.0;
	point2[2] = (point1[2] + point3[2])/2.0;
	point2[3] = (point1[3] + point3[3])/2.0;
	
	//Initialize values:
	A[1][1]=0.0, A[1][2]=0.0, A[1][3]=0.0;
	A[2][1]=0.0, A[2][2]=0.0, A[2][3]=0.0;
	A[3][1]=0.0, A[3][2]=0.0, A[3][3]=0.0;
	
	Z[1][1]=0.0;
	Z[2][1]=0.0;
	Z[3][1]=0.0;
	
	uknown[0]=0.0;
	uknown[1]=0.0;
	uknown[2]=0.0;
	//End INITIALIZATION.
	
	/*********************************************************************************
	 * Switching of the coordinates order because the parabola is x(y) = ay^2 + by + c
	 *********************************************************************************/
	/*if( !strcmp(function_type, "x(y)"))
	  {
	  //
	  //    ^x
	  //    |
	  //    |
	  //    |
	  //    |
	  //    |
	  //    ------------------> y
	
	    tmp = point1[1];
	    point1[1] = point1[2];
	    point1[2] = tmp;
	    tmp = point2[1];
	    point2[1] = point2[2];
	    point2[2] = tmp;
	    tmp = point3[1];
	    point3[1] = point3[2];
	    point3[2] = tmp;
	  }
	else if( !strcmp(function_type, "x(z)"))
	  {
	    //
	    //  ^x
	    //  |
	    //  |
	    //  |
	    //  |
	    //  |
	    //  ------------------> y
	    //
	    tmp = point1[1];
	    point1[1] = point1[3];
	    point1[3] = tmp;
	    tmp = point2[1];
	    point2[1] = point2[3];
	    point2[3] = tmp;
	    tmp = point3[1];
	    point3[1] = point3[3];
	    point3[3] = tmp;
	  }
	*/
	
	if( !strcmp(function_type, "z(x)"))
	  {
	    /*
	      ^z
	      |
	      |
	      |
	      |
	      |
	      ------------------> x
	    */
	    abscissa1 = point1[1];
	    ordinate1 = point1[3];
	    
	    abscissa2 = point2[1];
	    ordinate2 = point2[3];
	    
	    abscissa3 = point3[1];
	    ordinate3 = point3[3];

	  }
	else if( !strcmp(function_type, "y(x)"))
	  {
	    /*
	      ^y
	      |
	      |
	      |
	      |
	      |
	      ------------------> x
	    */
	    abscissa1 = point1[1];
	    ordinate1 = point1[2];
	    
	    abscissa2 = point2[1];
	    ordinate2 = point2[2];
	    
	    abscissa3 = point3[1];
	    ordinate3 = point3[2];

	  }

	if( !strcmp(function_type, "z(y)"))
	  {
	    /*
	      ^z
	      |
	      |
	      |
	      |
	      |
	      ------------------> y
	    */
	    abscissa1 = point1[2];
	    ordinate1 = point1[3];
	    
	    abscissa2 = point2[2];
	    ordinate2 = point2[3];
	    
	    abscissa3 = point3[2];
	    ordinate3 = point3[3];

	  }
	else if( !strcmp(function_type, "x(y)"))
	  {
	    /*
	      ^x
	      |
	      |
	      |
	      |
	      |
	      ------------------> y
	    */
	    abscissa1 = point1[2];
	    ordinate1 = point1[1];
	    
	    abscissa2 = point2[2];
	    ordinate2 = point2[1];
	    
	    abscissa3 = point3[2];
	    ordinate3 = point3[1];

	  }

	if( !strcmp(function_type, "x(z)"))
	  {
	    /*
	      ^x
	      |
	      |
	      |
	      |
	      |
	      ------------------> z
	    */
	    abscissa1 = point1[3];
	    ordinate1 = point1[1];
	    
	    abscissa2 = point2[3];
	    ordinate2 = point2[1];
	    
	    abscissa3 = point3[3];
	    ordinate3 = point3[1];

	  }
	else if( !strcmp(function_type, "y(z)"))
	  {
	    /*
	      ^y
	      |
	      |
	      |
	      |
	      |
	      ------------------> z
	    */
	    abscissa1 = point1[3];
	    ordinate1 = point1[2];
	    
	    abscissa2 = point2[3];
	    ordinate2 = point2[2];
	    
	    abscissa3 = point3[3];
	    ordinate3 = point3[2];

	  }

	
	/**************************************************
	 * Build the coefficient matrix and compute a,b,c:
	 **************************************************/
	A[1][1] = abscissa1*abscissa1;
	A[1][2] = abscissa1;
	A[1][3] = 1.0;
	
	A[2][1] = abscissa2*abscissa2;
	A[2][2] = abscissa2;
	A[2][3] = 1.0;
	
	A[3][1] = abscissa3*abscissa3;
	A[3][2] = abscissa3;
	A[3][3] = 1.0;
	
	Z[1][1] = ordinate1;
	Z[2][1] = ordinate2;
	Z[3][1] = ordinate3;
	
	dgaussj(A,3,Z,1);
				
	*a = Z[1][1];
	*b = Z[2][1];
	*c = Z[3][1];
	
	free_dmatrix(A,1,3,1,3);
	free_dmatrix(Z,1,3,1,1);
	free_dvector(point2,1,3);
	
return;
}

