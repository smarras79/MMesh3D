/************************************************************************************
 * ELLIPTIC solver for structured grids:
 *
 * Specifically, this is a Winslow (or smoothness) grid generator.
 * See [1]
 *
 * [1] P. Knupp and S. Steinberg, Fundamentals of grid generation. CRC Press
 *
 * simone.marras@gmail.com
 ************************************************************************************/

#include<stdio.h>
#include<math.h>

#include"myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

#define MAXITS_external 100
#define MAXITS 1000

#define EPS 1.0e-4

void elliptic_solver(double **COORDS, int **CONN, int *ELTYPE,		\
		     double **x, double **y, double **z, double **BDYFLAG, \
		     double **bottomSide, double **rightSide, double **topSide, double **leftSide, \
		     double *bottomSide_der, double *rightSide_der, double *topSide_der, double *leftSide_der, \
		     int nnodesx, int nnodesy, int nnodesz)
{
  
  int i,j,k;
  int iter;
  int nnodes, nelem;
  int internal_nodes;
 
  //double **BDYFLAG;
  double **a,**b,**c,**d,**e,**f;
  double *xbb, *ybb, *xbr, *ybr, *xbl, *ybl, *xbt, *ybt;
  double *Ubx, *Uby, *Utx, *Uty, *Vlx, *Vly, *Vrx, *Vry;
  double *XX, *YY;

  double xbb_der, ybb_der, xbr_der, ybr_der, xbl_der, ybl_der, xbt_der, ybt_der;
  double ybb1_der, ybbn_der, xbr1_der, xbrn_der, xbl1_der, xbln_der, ybt1_der, ybtn_der;
  double P, Q, apha, beta, H;
  
  //Hermite cubic polynomials:
  double H1psi,H2psi,H3psi,H4psi;
  double H1eta,H2eta,H3eta,H4eta;

  //Iteration values:
  double xij_prev, yij_prev, Fij, omega, err, tol;
  char *vtkfile_name;

  vtkfile_name = (char*) malloc(32 * sizeof(char *));

  nelem = (nnodesx - 1)*(nnodesy - 1);
  nnodes = nnodesx*nnodesy;
  
  f = dmatrix(1,nnodesx, 1,nnodesy);
    
  xbb = dvector(1,nnodesx), ybb = dvector(1,nnodesx);
  xbr = dvector(1,nnodesy), ybr = dvector(1,nnodesy);
  xbt = dvector(1,nnodesx), ybt = dvector(1,nnodesx);
  xbl = dvector(1,nnodesy), ybl = dvector(1,nnodesy);

  //Domain height:
  //H = ybt[1,1];
  //reminder:	bottomSide = dmatrix(1,nnodesx,1,2);

  //Build connectivity matrix:
  GRID2CONN(0, 0, nnodesx, nnodesy, CONN, ELTYPE, 1);

  /*
   * Global cycle:
   */
  k = 1;
  tol = 1e-4;
  
  iter = 0;
  for( iter = 1; iter<MAXITS_external; iter++)
    {
      
      ///B.C.
      for(i=1; i<=nnodesx; i++){
	x[i][1] = bottomSide[i][1];
	y[i][1] = bottomSide[i][2];
	
	x[i][nnodesy] = topSide[i][1];
	y[i][nnodesy] = topSide[i][2]; 
      } 
      
      for(j=1; j<=nnodesy; j++){
	x[nnodesx][j] = rightSide[j][1];
	y[nnodesx][j] = rightSide[j][2];
	
	x[1][j] = leftSide[j][1];
	y[1][j] = leftSide[j][2];
      }
      
      //Solve the system with Gauss-Seidel:
      gauss_seidel_xy( x, y, z, nnodesx, nnodesy, nnodesz, &err);
      
      //Store VTK at every iteration:
      //Store COORDS:
      
      k = 0;
      for (j=1; j<=nnodesy; j++)
	for (i=1; i<=nnodesx; i++)
	  {
	    k = k + 1;
	    
	    //Store Coords:
	    COORDS[k][1] = k;
	    COORDS[k][2] = x[i][j];
	    COORDS[k][3] = y[i][j];
	  }
      
      //Print to VTK at every iteration:
      sprintf(vtkfile_name, "mesh_vtk%d.vtk", iter);
      dwrt2VTK(vtkfile_name);
      
      if( err <= tol )
	{
	  //printf(" Converged at Iteration: %d\n", iter);
	  continue;
	}

    }
  
  //Store COORDS:
  k = 0;
  for (j=1; j<=nnodesy; j++)
    for (i=1; i<=nnodesx; i++)
      {
	k = k + 1;

	  //Store Coords:
	COORDS[k][1] = k;
	COORDS[k][2] = x[i][j];
	COORDS[k][3] = y[i][j];
	
	//Assign the boundary flags:
	if( j == 1 && i>1 && i<nnodesx){
	  
	  BDYFLAG[k][1] = 1;				// == 1 means that is it a boundary node on the BOTTOM boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( i == nnodesx && j>1 && j<nnodesy){
	  
	  BDYFLAG[k][1] = 2;				// == 2 means that is it a boundary node on the EAST boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( j == nnodesy && i>1 && i<nnodesx){
	  
	  BDYFLAG[k][1] = 3;				// == 3 means that is it a boundary node on the TOP boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( i == 1 && j>1 && j<nnodesy){
	  
	  BDYFLAG[k][1] = 4;				// == 4 means that is it a boundary node on the WEST boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( j == 1 && i == nnodesx){
	  
	  BDYFLAG[k][1] = 5;				// == 5 means that is it a boundary node on the BOTTOM-RIGHT boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( i == nnodesx && j == nnodesy){
	  
	  BDYFLAG[k][1] = 6;				// == 6 means that is it a boundary node on the TOP-RIGHT boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( j == nnodesy && i == 1){
	  
	  BDYFLAG[k][1] = 7;				// == 7 means that is it a boundary node on the TOP-LEFT boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else if( i == 1 && j == 1){
	  
	  BDYFLAG[k][1] = 8;				// == 8 means that is it a boundary node on the BOTTOM-LEFT boundary.
	  BDYFLAG[k][2] = COORDS[k][2];
	  BDYFLAG[k][3] = COORDS[k][3];
	}
	else{
	  
	  BDYFLAG[k][1] = 0;				// == 0 means that is it NOT a boundary node.
	  BDYFLAG[k][2] = 99999999999.0;
	  BDYFLAG[k][3] = 99999999999.0;
	}
	
	// Store the boundary nodes:
      }
  
  //Free memory:
  free(vtkfile_name);
  free_dmatrix(f, 1,nnodesx, 1, nnodesy);
   
  free_dvector(xbb, 1,nnodesx), free_dvector(ybb, 1,nnodesx);
  free_dvector(xbr, 1,nnodesy), free_dvector(ybr, 1,nnodesy);
  free_dvector(xbt, 1,nnodesx), free_dvector(ybt, 1,nnodesx);
  free_dvector(xbl, 1,nnodesy), free_dvector(ybl, 1,nnodesy);
  
  return; //BDYFLAG;
}//END elliptic_solver

void gauss_seidel_xy( double **x, double **y, double **z, int nnodesx, int nnodesy, int nnodesz, double *err)
{
  /*
   * 1 iteration of Gauss-Seidel:
   */
  int i, j;
  double xtemp, ytemp;
  double g11, g22, g12;
  double dksi,dksi2, deta,deta2;

  dksi = 1.0/(nnodesx-1);
  deta = 1.0/(nnodesy-1);
 
  dksi2 = dksi*dksi;
  deta2 = deta*deta;

  *err = 0.0;
  for( j = 2; j < nnodesy; j++)   
    for( i = 2; i < nnodesx; i++)
      {
	g11 = ( (x[i][j+1] - x[i][j-1])*(x[i][j+1] - x[i][j-1])/dksi2 +
		(y[i][j+1] - y[i][j-1])*(y[i][j+1] - y[i][j-1])/deta2 )/4;
	
	g22 = ( (x[i+1][j] - x[i-1][j])*(x[i+1][j] - x[i-1][j])/dksi2 +
		(y[i+1][j] - y[i-1][j])*(y[i+1][j] - y[i-1][j])/deta2 )/4;
	
	g12 = ( (x[i+1][j] - x[i-1][j])*(x[i][j+1] - x[i][j-1])/dksi2 + 
		(y[i+1][j] - y[i-1][j])*(y[i][j+1] - y[i][j-1])/deta2 )/(4*dksi*deta);
	
	xtemp = dksi2/(2*(g11+g22)) * (g11*(x[i+1][j] + x[i-1][j])/dksi2 + 
				       g22*(x[i][j+1] + x[i][j-1])/dksi2
				       - 0.5*g12*x[i+1][j+1] + 0.5*g12*x[i-1][j+1] +
				       - 0.5*g12*x[i-1][j-1] + 0.5*g12*x[i+1][j-1]);
	
	ytemp = deta2/(2*(g11+g22)) * (g11*(y[i+1][j] + y[i-1][j])/deta2 + 
				       g22*(y[i][j+1] + y[i][j-1])/deta2
				       - 0.5*g12*y[i+1][j+1] + 0.5*g12*y[i-1][j+1] +
				       - 0.5*g12*y[i-1][j-1] + 0.5*g12*y[i+1][j-1]);
      
	*err += (x[i][j]-xtemp)*(x[i][j]-xtemp) + (y[i][j]-ytemp)*(y[i][j]-ytemp);
	//printf("  err=%lf\n",  *err);
	x[i][j] = xtemp;
	y[i][j] = ytemp;

      }
  //printf(" ind =%d\n", ind);
  *err = sqrt( *err/((nnodesx-2)*(nnodesy-2)) );
  
  return;
}//End GS

/*
 * Compute f (distorsion function) as in Eca, (1996)
 */

void compute_f( double **x, double **y, int nnodesx, int nnodesy, double **f)
{
  int i,j;
  double a,b,c,d;
  
  for( i=2; i<nnodesx; i++)
      for( j=2; j<nnodesy; j++)
	{
	  //f(i+1/2,j)
	  a = (x[i+1][j+1] + x[i][j+1] - x[i+1][j-1] - x[i][j-1])*(x[i+1][j+1] + x[i][j+1] - x[i+1][j-1] - x[i][j-1]);
	  b = (y[i+1][j+1] + y[i][j+1] - y[i+1][j-1] - y[i][j-1])*(y[i+1][j+1] + y[i][j+1] - y[i+1][j-1] - y[i][j-1]);
	  c = (x[i+1][j] - x[i][j])*(x[i+1][j] - x[i][j]);
	  d = (y[i+1][j] - y[i][j])*(y[i+1][j] - y[i][j]);
 
	  f[i+1][j] = 0.5*sqrtf( (a + b)/(c + d) );


	  //f(i-1/2,j)
	  a = (x[i-1][j+1] + x[i][j+1] - x[i-1][j-1] - x[i][j-1])*(x[i-1][j+1] + x[i][j+1] - x[i-1][j-1] - x[i][j-1]);
	  b = (y[i-1][j+1] + y[i][j+1] - y[i-1][j-1] - y[i][j-1])*(y[i-1][j+1] + y[i][j+1] - y[i-1][j-1] - y[i][j-1]);
	  c = (x[i][j] - x[i-1][j])*(x[i][j] - x[i-1][j]);
	  d = (y[i][j] - y[i-1][j])*(y[i][j] - y[i-1][j]);

	  f[i-1][j] = 0.5*sqrtf( (a + b)/(c + d) );


	  //f(i,j+1/2)
	  a = (x[i+1][j] + x[i+1][j+1] - x[i-1][j] - x[i-1][j+1])*(x[i+1][j] + x[i+1][j+1] - x[i-1][j] - x[i-1][j+1]);
	  b = (y[i+1][j] + y[i+1][j+1] - y[i-1][j] - y[i-1][j+1])*(y[i+1][j] + y[i+1][j+1] - y[i-1][j] - y[i-1][j+1]);
	  c = (x[i][j+1] - x[i][j])*(x[i][j+1] - x[i][j]);
	  d = (y[i][j+1] - y[i][j])*(y[i][j+1] - y[i][j]);

	  f[i][j+1] = 0.5*sqrtf( (c + d)/(a + b) );


	  //f(i,j-1/2)
	  a = (x[i+1][j] + x[i+1][j-1] - x[i-1][j] - x[i-1][j+1])*(x[i+1][j] + x[i+1][j-1] - x[i-1][j] - x[i-1][j+1]);
	  b = (y[i+1][j] + y[i+1][j-1] - y[i-1][j] - y[i-1][j+1])*(y[i+1][j] + y[i+1][j-1] - y[i-1][j] - y[i-1][j+1]); 
	  c = (x[i][j] - x[i][j-1])*(x[i][j] - x[i][j-1]);
	  d = (y[i][j] - y[i][j-1])*(y[i][j] - y[i][j-1]);

	  f[i][j-1] = 0.5*sqrtf( (c + d)/(a + b) );

	}

  return;
}

/*
 * Iterative solver of the Thompson-Thames-Mastin elliptic system:
 *
 * Simone Marras, January 2012
 */
void apply_smoothing()
{
  int i, j, k;
  int ii,jj,kk;
  int iter;
  double err;
  
  MEMORY_ALLOCATE(4);

  ii = 1;
  for(iter=1; iter<=20; iter++)
    {
      for(k=1; k<=nnodesz; k++)
	for(j=1; j<=nnodesy; j++)
	  for(i=1; i<=nnodesx; i++)
	    {
	      x[i][j][k] = COORDS[ii][2];
	      y[i][j][k] = COORDS[ii][3];
	      z[i][j][k] = COORDS[ii][4];
	      
	      ii = ii + 1;
	    }
      //      gauss_seidel_TTM(iter, x, y, z, err);
    }
  
  //MEMORY_DEALLOCATE(4);
}

void gauss_seidel_TTM(int iter, double ***x, double ***y, double ***z, double err)
{ 
  /*
   * 1 iteration of Gauss-Seidel:
   */
  
  int i, j, k, ipoin2d, ipoin3d;
  double  xtemp, ytemp, ztemp, rhs;
  double  g11, g22, g33, g12, g13, g23, Jac,Jinv, J;

  double  P, Q, R;
  double  tol;
  double  dx, dy, dz;

  double  e_x,e_y,e_z;
  double  n_x,n_y,n_z;
  double  c_x,c_y,c_z;

  double  x_e, x_n, x_c;
  double  y_e, y_n, y_c;
  double  z_e, z_n, z_c;
  
  double  x_en, x_ec, x_nc;
  double  y_en, y_ec, y_nc;
  double  z_en, z_ec, z_nc;
  
  /*
   * If not defined, set the lowest layer to 1 
   * (first iz in the vertical direction)
   */
  //if(lowest_layer < 1) lowest_layer = 1;
  
  /*
   * Solve the Thompso--Thames-Mastin system using Gauss-Seidel
   */
  //for(k=lowest_layer+1; k<=nnodesz-1; k++)
  for(k=1; k<=nnodesz-1; k++)
    for(j=2; j<=nnodesy-1; j++)
      for(i=2; i<=nnodesx-1; i++)
	{
	  //First derivatives:
	  x_e = 0.5*(x[i+1][j][k] - x[i-1][j][k]);
	  x_n = 0.5*(x[i][j+1][k] - x[i][j-1][k]);
	  x_c = 0.5*(x[i][j][k+1] - x[i][j][k-1]);
  
	  y_e = 0.5*(y[i+1][j][k] - y[i-1][j][k]);
	  y_n = 0.5*(y[i][j+1][k] - y[i][j-1][k]);
	  y_c = 0.5*(y[i][j][k+1] - y[i][j][k-1]);

	  z_e = 0.5*(z[i+1][j][k] - z[i-1][j][k]);
	  z_n = 0.5*(z[i][j+1][k] - z[i][j-1][k]);
	  z_c = 0.5*(z[i][j][k+1] - z[i][j][k-1]);

	  //Crossed second derivatives:
	  x_en = 0.25*( x[i+1][j+1][k] - x[i-1][j+1][k] + x[i-1][j-1][k] - x[i+1][j-1][k] );
	  x_ec = 0.25*( x[i+1][j][k+1] - x[i-1][j][k+1] + x[i-1][j][k-1] - x[i+1][j][k-1] );
	  x_nc = 0.25*( x[i][j+1][k+1] - x[i][j-1][k+1] + x[i][j-1][k-1] - x[i][j+1][k-1] );

	  y_en = 0.25*( y[i+1][j+1][k] - y[i-1][j+1][k] + y[i-1][j-1][k] - y[i+1][j-1][k] );
	  y_ec = 0.25*( y[i+1][j][k+1] - y[i-1][j][k+1] + y[i-1][j][k-1] - y[i+1][j][k-1] );
	  y_nc = 0.25*( y[i][j+1][k+1] - y[i][j-1][k+1] + y[i][j-1][k-1] - y[i][j+1][k-1] );

	  z_en = 0.25*( z[i+1][j+1][k] - z[i-1][j+1][k] + z[i-1][j-1][k] - z[i+1][j-1][k] );
	  z_ec = 0.25*( z[i+1][j][k+1] - z[i-1][j][k+1] + z[i-1][j][k-1] - z[i+1][j][k-1] );
	  z_nc = 0.25*( z[i][j+1][k+1] - z[i][j-1][k+1] + z[i][j-1][k-1] - z[i][j+1][k-1] );

	  //Jacobian term:
	  Jac = x_e*y_n*z_c + x_n*y_c*z_e + x_c*y_e*z_n - x_e*y_c*z_n - x_n*y_e*z_c - x_c*y_n*z_e;
	  Jinv = 1.0/J;
           
	  //Back derivatives:
	  e_x = Jinv*(y_n*z_c - y_c*z_n);
	  e_y = Jinv*(x_c*z_n - x_n*z_c);
	  e_z = Jinv*(x_n*y_c - x_c*y_n);

	  n_x = Jinv*(y_c*z_e - y_e*z_c);
	  n_y = Jinv*(x_e*z_c - x_c*z_e);
	  n_z = Jinv*(x_c*y_e - x_e*y_c);

	  c_x = Jinv*(y_e*z_n - y_n*z_e);
	  c_y = Jinv*(x_n*z_e - x_e*z_n);
	  c_z = Jinv*(x_e*y_n - x_n*y_e);

	  //Metric tensor G (components of -):
	  g11 = e_x*e_x + e_y*e_y + e_z*e_z;
	  g22 = n_x*n_x + n_y*n_y + n_z*n_z;
	  g33 = c_x*c_x + n_y*c_y + c_z*c_z;
	  g12 = e_x*n_x + e_y*n_y + e_z*n_z;
	  g13 = e_x*c_x + e_y*c_y + e_z*c_z;
	  g23 = n_x*c_x + n_y*c_y + n_z*c_z;
	     
	  //Compute the control functions P and Q:
	  //	  compute_PQ( i, j, k, dx, dy, dz, P, Q, R);
	  
	  //X:
	  rhs =						\
	    g11*(x[i-1][j][k] + x[i+1][j][k]) +		\
	    g22*(x[i][j-1][k] + x[i][j+1][k]) +		\
	    g33*(x[i][j][k-1] + x[i][j][k+1]) +		\
	    2.0*(g12*x_en + g13*x_ec + g23*x_nc);
	  
	  xtemp = 0.5*rhs/(g11 + g22 + g33);
	  
	  //Y:
	  rhs =						\
	    g11*(y[i-1][j][k] + y[i+1][j][k]) +		\
	    g22*(y[i][j-1][k] + y[i][j+1][k]) +		\
	    g33*(y[i][j][k-1] + y[i][j][k+1]) +		\
	    2.0*(g12*y_en + g13*y_ec + g23*y_nc);
	  
	  ytemp = 0.5*rhs/(g11 + g22 + g33);
	  
	  //Z:
	  rhs =					   \
	    g11*(z[i-1][j][k] + z[i+1][j][k]) +	   \
	    g22*(z[i][j-1][k] + z[i][j+1][k]) +	   \
	    g33*(z[i][j][k-1] + z[i][j][k+1]) +	   \
	    2.0*(g12*z_en + g13*z_ec + g23*z_nc) - \
	    g33*z_c*R;
	  
	  ztemp = 0.5*rhs/(g11 + g22 + g33);
	  
	  //Adjust error at current iteration:
	  err = err +					\
	    (x[i][j][k]-xtemp)*(x[i][j][k]-xtemp) +	\
	    (y[i][j][k]-ytemp)*(y[i][j][k]-ytemp) +	\
	    (z[i][j][k]-ztemp)*(z[i][j][k]-ztemp);
	  
	  x[i][j][k] = xtemp;
	  y[i][j][k] = ytemp;
	  z[i][j][k] = ztemp;
	}
  
    err = sqrt( err/((nnodesx-2)*(nnodesy-2)*(nnodesz-2)));
  
    return;
}
