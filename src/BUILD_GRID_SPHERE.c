/*
 * This function builds the volume grid on a spherical surface
 *
 * The original spherical HEXA and ICO grids were written by F. X. Giraldo within NUMA3DCG.
 *
 * C adaptation of the original F90 code translated by:
 * Simone Marras, April 2014
 *
 */
#include "myinclude.h"
#include "mydefine.h"
#include "global_vars.h"

int BUILD_GRID_SPHERE(int rank)
{
  int nel;
  int nelem_r, nelem_s;
  int npoin_r, npoin_s;
  
  int *intma_r;
  int *intma_s;
  
  double radius;
  double rmin, rmax;

  nel     = nelx;
  nelem_r = nelz;
  radius  = earth_radius;
  
  if( !strcmp(problem[0], "hexa") || !strcmp(problem[0], "hex") )
    {
      nelem_s = 6*nel*nel;                 /* Number of elements on the sphere surface */
      nelem_g = nelem_s*nelem_r;           /* Total number of elements */
      npoin_r = nopz*nelem_r + 1;          /* Number of points in the radial direction */
      npoin_s = 6*(nel*nop)*(nel*nop) + 2; /* Number of points on the sphere surface. We are taking nopy=nopx=nop */

      if(rank == 0) 
	{
	  printf(" nelem_s = %d\n", nelem_s);
	  printf(" nelem_g = %d\n", nelem_g);
	  printf(" npoin_r = %d\n", npoin_r);
	  printf(" npoin_s = %d\n", npoin_s);

	}
      
      /* Allocate local memory */
      intma_r = ivector(1,nglz*nelem_r);
      intma_s = ivector(1,nglx*ngly*nelem_s);
      
      /* Create grid: */
      create_grid_sphere_hex(COORDS, CONN, intma_r,		\
			     intma_s,				\
			     nel, nelem_r, nelem_s,		\
			     npoin_r, npoin_s, rmin, rmax);
	
	
      /* Free local memory */
      free_ivector(intma_r, 1,nglz*nelem_r);
      free_ivector(intma_s, 1,nglx*ngly*nelem_s);

    }
  
  else if ( !strcmp(problem[0], "icosahedron") || !strcmp(problem[0], "ico") )
    {
      nelem_s = 60*nel*nel;
      nelem_g =  nelem_s*nelem_r; 
      npoin_r = nopz*nelem_r + 1;     /* Number of points in radial direction */
      npoin_s = 60*(nel*nopx)*(nel*nopx) + 2; /* Number of points on each sphere: nopy=nopx */

      if(rank == 0) 
	{
	  printf(" nelem_s = %d\n", nelem_s);
	  printf(" nelem_g = %d\n", nelem_g);
	  printf(" npoin_r = %d\n", npoin_r);
	  printf(" npoin_s = %d\n", npoin_s);
	}
      
    }
  
  nnodes_g = npoin_s*npoin_r;   /* Total number of points */
  nboun_g = 2*nelem_s;         /* Number of bounding elements */
  ncol_g  = nnodes_g/npoin_r;
  
  return 0;
}

int create_grid_sphere_hex(double **coord, int **intma, int *intma_r, \
			   int ***intma_s, int nel, int nelem_r, int nelem_s, \
			   int npoin_r, int npoin_s, double rmin, double rmax)
/*int create_grid_sphere_hex(double **coord, int **intma, int *intma_r, \
  int *intma_s, int *ele_col, int **bsido,	\
  int *iboun, int nel, int nelem_r, int nelem_s, \
  int npoin_r, int npoin_s, double rmin, double rmax)
*/
{
  /*
    double coord(3,npoin)
    int indexg(2,npoin)
    int intma(nglx,ngly,nglz,nelem)
    int intma_lev(nglx,ngly,nelem_s,npoin_r)
    int ele_col(nelem)
    int bsido(6,nboun), iboun(2)
    int nboun, npoin_r, npoin_s,nel
    int nglx,ngly,nglz, nelem,nelem_r,nelem_s,npoin
    double rmin, rmax
  */
  int ier, ipr, ies, ips, ie, ip, ib;
  int i, j, k;
  double dr;
  double r2;
  
  double  **coord_s;
  double   *r;
  
  int nface;
  int nx, ny;
  int npoin0;
  int nelem0;

  nface  = 6;
  nx     = nel*(nglx-1) + 1;
  ny     = nel*(ngly-1) + 1;
  npoin0 = nx*ny;
  nelem0 = nel*nel;

  rmin = earth_radius;
  rmax = rmin + zmax;
  
  /* Allocate local memory and initialize to zero */
  r       = dvector(1,npoin_r);
  for (i=1; i<=npoin_r; i++)
    r[i] = 0.0;
  
  coord_s = dmatrix(1,3, 1,nnodes);
  for (i=1; i<=nnodes; i++)
    {
      coord_s[1][i] = 0.0;
      coord_s[2][i] = 0.0;
      coord_s[3][i] = 0.0;
    }
  
  
  /* Create 1D Spectral Element Grid in the Radial direction */
  create1d_grid(r,intma_r, nelem_r, npoin_r, rmin, rmax);
  
  /* Create 2D Spherical Grid: assumes nely=nelx */
  hex_quad(coord_s,intma_s,nel,nglx,nface,nelem_s,npoin_s,nelem0,npoin0,nx,ny);
      
/*
    !Construct Cartesian Product [rmin,rmax] x S^2 to create
    !a spherical shell
    do ier = 1,nelem_r
    do k = 1,nglz
        
    ipr = intma_r(k,ier)
    do ies = 1,nelem_s
    do i = 1,nglx
    do j = 1,ngly
    ips = intma_s(i,j,ies)
    ie = (ier - 1)*nelem_s + ies
    ip = (ipr - 1)*npoin_s + ips
    intma(i,j,k,ie) = ip                !3D Intma
    intma_lev(i,j,ies,ipr) = ip         !Level-Based Intma
    ele_col(ie) = ies                   !Level-Based element that belongs to ie
    coord(1,ip) = r(ipr)*coord_s(1,ips)
    coord(2,ip) = r(ipr)*coord_s(2,ips)
    coord(3,ip) = r(ipr)*coord_s(3,ips)
    indexg(1,ip) = ips
    indexg(2,ip) = ipr
    end do
    end do
    end do
    end do
    end do
    print *, "Coords Created"
  
    !Create Bsido
    !NFBC on the Ground
    ier = 1
    ib = 1
    do ies =1,nelem_s
    ie = (ier - 1)*nelem_s + ies 
    k = 1
     
    i = 1; j = 1
    ip = intma(i,j,k,ie)
    bsido(1,ib) = ip
     
    i = 1; j = ngly
    ip = intma(i,j,k,ie)
    bsido(2,ib) = ip
     
    i = nglx; j = ngly
    ip = intma(i,j,k,ie)
    bsido(3,ib) = ip
     
    i = nglx; j = 1
    ip = intma(i,j,k,ie)
    bsido(4,ib) = ip
    bsido(5,ib) = ie
    bsido(6,ib) = iboun(1)
     
    ib = ib + 1
    end do
  
    !NRBC at the top of the atmosphere        
    ier = nelem_r
    do ies =1,nelem_s
    ie = (ier - 1)*nelem_s + ies
    k = nglz
     
    i = 1; j = 1
    ip = intma(i,j,k,ie)
    bsido(1,ib) = ip
     
    i = nglx; j = 1
    ip = intma(i,j,k,ie)
    bsido(2,ib) = ip
     
    i = nglx; j = ngly
    ip = intma(i,j,k,ie)
    bsido(3,ib) = ip
     
    i = 1; j = ngly
    ip = intma(i,j,k,ie)
    bsido(4,ib) = ip
    bsido(5,ib) = ie
    bsido(6,ib) = iboun(2)
     
    ib = ib + 1
     
    end do
  */


  /* Deallocate local memory */
  free_dvector(r, 1,npoin_r);
  free_dmatrix(coord_s, 1,3, 1,nnodes);
  
  return 0;
}

/*
 * Create 1D grid:
 */
int create1d_grid(double *r, int *intmar, int nelr, int npoinr, double rmin, double rmax)
{
  
  int i;
  int ier, ipr;
  
  double num, den;
  
  double *z;
  double *dr;
  
  /* Allocate local arrays */
  z  = dvector(1,nelr+1);
  dr = dvector(1,nelr+1);
  
  /* Initialize to zero */
  for(i=1; i<=nelr+1; i++)
    {
      z[i]  = 0.0;
      dr[i] = 0.0;
    }
  
  for(i=2; i <= nelr + 1; i++)
    dr[i-1] = (rmax - rmin)/nelr;
  
  /* Compute the LGL points */
  legendre_gauss_lobatto(ngl,xgl,wgl);
  
  ier = 1;
  ipr = 0;
  for(i=1; i<=ngl; i++)
    {
      ipr = ipr + 1;
      
      r[ipr] = rmin +  0.5*dr[ier]*(xgl[i] + 1.0);
      intmar[i*ier] = ipr;
    }
  
  z[1] = rmin;
  for(ier=2; ier<=nelr; ier++)
    {
      z[ier] = z[ier-1] + dr[ier];
      
      //intmar(1,ier) = ipr
      intmar[ier] = ipr;
      for( i=2; i<=ngl; i++)
	{
	  ipr = ipr + 1;
	  r[ipr] = z[ier] + dr[ier]*(xgl[i] + 1.0)*0.5;

	  //intmar(i,ier) = ipr;
	  intmar[i*ier] = ipr;
	  
	}
    }
  
  /* Free local memory */
  free_dvector(z,  1,nelr+1);
  free_dvector(dr, 1,nelr+1);
  
  return 0;
}


/*----------------------------------------------------------------------!
  !This subroutine constructs the HEXAHEDRAL GRID
  !Written by Francis X. Giraldo on 10/98
  !           Naval Research Laboratory
  !           Global Modeling Section
  !           Monterey, CA 93943-5502
  !----------------------------------------------------------------------!
*/
int hex_quad(double **coord, int **intma, int nel, int ngl, int nface, int nelem, \
	     int npoin, int nelem0, int npoin0, int nx, int ny)
{  
  double  **xcoors0,  **ycoors0;
  double  **xcoord0,  **ycoord0;
  double ***xcoorsf, ***ycoorsf;
  double ***xcoordf, ***ycoordf, ***zcoordf;
  double   *xcoors,    *ycoors,    *zcoors;
  
  int *intma0;
  int *nodef;
  
  /* Face Information */
  double *lon, *lat;
  
  double pi, twopi, pio2;
  double clon, clat, alon, alat, rlon, rlat, r;
  double x, y, z, b_gnom1, b_gnom2, b_rot1, b_rot2;
  int i, j, k, l, m, ip, np, ie, ii, jj;
  
  /* Allocate local memory */
  xcoors0 = dmatrix(1,nx,1,ny); ycoors0 = dmatrix(1,nx,1,ny);
  xcoord0 = dmatrix(1,nx,1,ny); ycoord0 = dmatrix(1,nx,1,ny);
  
  xcoorsf = d3tensor(1,nx,1,ny,1,nface); ycoorsf = d3tensor(1,nx,1,ny,1,nface);
  xcoordf = d3tensor(1,nx,1,ny,1,nface); ycoordf = d3tensor(1,nx,1,ny,1,nface); zcoordf = d3tensor(1,nx,1,ny,1,nface);
  xcoors = dvector(1,npoin), ycoors = dvector(1,npoin), zcoors = dvector(1,npoin);
  
  lon = dvector(1,6);
  lat = dvector(1,6);

  lon[1] = 0.0;
  lon[2] = 90.0;
  lon[3] = 180.0;
  lon[4] = 270.0;
  lon[5] = 0.0;
  lon[6] = 0.0;
  
  lat[1] = 0.0;
  lat[2] = 0.0;
  lat[3] = 0.0;
  lat[5] = 0.0;
  lat[6] = 90.0;
  lat[7] = -90.0;
  
  intma0 = ivector(1,ngl*ngl*nelem0);
  nodef = ivector(1,nx*ny*nface);
  
  /* Constants */
  pi    = 4.0*atan(1.0);
  twopi = 2*pi;
  pio2  = 0.5*pi;
  
  for (i=1; i<=6; i++)
    {
      lon[i] = lon[i]*pi/180.0;
      lat[i] = lat[i]*pi/180.0;
    }
  
  /* Construct the 2D Planar Grid */
  init_hex(xcoord0,ycoord0,intma0,nelem0,nel,nx,ny,ngl);
  
  /*
    !Get Rotated Spherical Coordinates
    do j=1,ny
    do i=1,nx
    x=coord0(1,i,j)
    y=coord0(2,i,j)
        
    !Get Backward Gnomonic Projection
    rlon=b_gnom1(x)
    rlat=b_gnom2(y,rlon)
        
    !Store Spherical coordinates
    coors0(1,i,j)=rlon
    coors0(2,i,j)=rlat
        
    end do !i   
    end do !j   
  
    !Loop through the 6 Faces
    ip=0
    do l=1,nface
    clon=lon(l)
    clat=lat(l)
     
    !loop through the Planar Grid and get the Face Coordinates
    do j=1,ny
    do i=1,nx
           
    !Get Rotated Spherical Coordinates
    rlon=coors0(1,i,j)
    rlat=coors0(2,i,j)
           
    !Get Backward Rotation Transformation
    alon=b_rot1(rlon,rlat,clon,clat)
    alat=b_rot2(rlon,rlat,clon,clat)
           
    if (alon < 0 ) alon=alon + twopi
    if (alon > twopi ) alon=alon - twopi
           
    !Store Spherical coordinates
    coorsf(1,i,j,l)=alon
    coorsf(2,i,j,l)=alat
           
    end do !i   
    end do !j   
    end do !l
  
    !!--------Unite Faces into Global Grid---------!!
  
    np=0
  
    !FACE 1
    do j=1,ny
    do i=1,nx
    np=np + 1
    nodef(i,j,1)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,1)
    end do !k
    end do !i
    end do !j   
  
    !FACE 2
    do j=1,ny
    nodef(1,j,2)=nodef(nx,j,1)
    end do !j   
    do j=1,ny
    do i=2,nx
    np=np + 1
    nodef(i,j,2)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,2)
    end do !k
    end do !i
    end do !j   
  
    !FACE 3
    do j=1,ny
    nodef(1,j,3)=nodef(nx,j,2)
    end do !j   
    do j=1,ny
    do i=2,nx
    np=np + 1
    nodef(i,j,3)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,3)
    end do !k
    end do !i
    end do !j   
  
    !FACE 4
    do j=1,ny
    nodef(1,j,4)=nodef(nx,j,3)
    end do !j   
    do j=1,ny
    nodef(nx,j,4)=nodef(1,j,1)
    end do !j   
    do j=1,ny
    do i=2,nx-1
    np=np + 1
    nodef(i,j,4)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,4)
    end do !k
    end do !i
    end do !j   
  
    !FACE 5
    do i=1,nx
    nodef(i,1,5)=nodef(i,ny,1)
    nodef(i,ny,5)=nodef(nx+1-i,ny,3)
    end do !j   
    do j=1,ny
    nodef(1,j,5)=nodef(nx+1-j,ny,4)
    nodef(nx,j,5)=nodef(j,ny,2)
    end do !j   
    do j=2,ny-1
    do i=2,nx-1
    np=np + 1
    nodef(i,j,5)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,5)
    end do !k
    end do !i
    end do !j   
  
    !FACE 6
    do i=1,nx
    nodef(i,1,6)=nodef(nx+1-i,1,3)
    nodef(i,ny,6)=nodef(i,1,1)
    end do !j   
    do j=1,ny
    nodef(1,j,6)=nodef(j,1,4)
    nodef(nx,j,6)=nodef(nx+1-j,1,2)
    end do !j   
    do j=2,ny-1
    do i=2,nx-1
    np=np + 1
    nodef(i,j,6)=np
    do k=1,2
    coors(k,np)=coorsf(k,i,j,6)
    end do !k
    end do !i
    end do !j   
  
    if (np /= npoin) then 
    print*,' npoin np = ',npoin,np
    stop
    end if
  
    !GENERATE INTMA
    ie=0
    do m=1,nface
    do l=1,nel
    do k=1,nel
    ie=ie+1
    do j=1,ngl
    jj=(ngl-1)*(l-1) + j
    do i=1,ngl
    ii=(ngl-1)*(k-1) + i
    ip=nodef(ii,jj,m)
    intma(i,j,ie)=ip
    end do !i
    end do !j
    end do !k   
    end do !l
    end do !m
  
    !Convert from Spherical to Cartesian
    do i=1,npoin
    alon=coors(1,i)
    alat=coors(2,i)
    x=cos(alat)*cos(alon)
    y=cos(alat)*sin(alon)
    z=sin(alat)
    r=sqrt(x*x + y*y + z*z)
    coord(1,i)=x/r
    coord(2,i)=y/r
    coord(3,i)=z/r
     
    end do !i   
  */

  /* Free local memory */
  free_dvector(lon,1,6);
  free_dvector(lat,1,6);
  
  free_ivector(intma0, 1,ngl*ngl*nelem0);
  free_ivector(nodef, 1,nx*ny*nface);
  
  return 0;
}//end function hex_quad

/*!----------------------------------------------------------------------!
  !This subroutine constructs the INITIAL HEXAHEDRON
  !Written by Francis X. Giraldo on 10/98
  !           Naval Research Laboratory
  !           Global Modeling Section
  !           Monterey, CA 93943-5502
  !----------------------------------------------------------------------!
*/
int init_hex(double **xcoord, double **ycoord, int **intma, int nelem, int nel, int nx, int ny, int ngl)
{
  //local arrays
  double *xgl, *wgl;

  int **node;
  
  double pi;
  double xmin, xmax;
  double ymin, ymax;
  double dx, dy;
  double xl, yl;
  double x0, y0;
  double x, y;

  int i, k;
  int ip, ii, jj;
  int l, l1, l2;
  int j, j1, j2;

  //Allocate local arrays
  xgl = dvector(1,ngl);
  wgl = dvector(1,ngl);
  node = imatrix(1,nx, 1,ny);

  //set some constants
  pi=PI;
  xmin=-1.0;
  xmax=+1.0;
  ymin=-1.0;
  ymax=+1.0;
  dx=(xmax-xmin)/(nel);
  dy=(ymax-ymin)/(nel);
  xl=xmax-xmin;
  yl=ymax-ymin;
  
  //Generate Gauss-Lobatto Points for NPTS
  legendre_gauss_lobatto(ngl,xgl,wgl);
  
  //GENERATE COORD
  ip=0;
  jj=0;
  for( k=1; k<=nel; k++)
    {
      y0=ymin + (double)(k-1)*dy;
     
      if (k == 1)
	l1=1;
      else
	l1=2;
      
      for( l=l1; l<=ngl; l++)
	{
	  y=( xgl[l] + 1.0 )*dy*0.5 + y0;
	 
	  jj=jj+1;
	  ii=0;
	 
	  for( i=1; i<=nel; i++)
	    {
	      x0=xmin + (double)(i-1)*dx;
	     
	      if (i == 1) 
		j1=1;
	      else
		j1=2;
	      
	      for( j=j1; j<=ngl; j++)
		{
		  ii=ii + 1;
		  ip=ip + 1;
		  x=( xgl[j]+1 )*dx*0.5 + x0;
		  xcoord[ii][jj] = x;
		  ycoord[ii][jj] = y;
		  node[ii][jj]   = ip;
		} 
	    } 
	}
    }
  
  /*
    !GENERATE INTMA
    ie=0
    do k=1,nel
    do i=1,nel
    ie=ie+1
    do l=1,ngl
    jj=(ngl-1)*(k-1) + l
    do j=1,ngl
    ii=(ngl-1)*(i-1) + j
    ip=node(ii,jj)
    intma(j,l,ie)=ip
    end do
    end do
    end do
    end do
  */


  //Deallocate local arrays
  free_dvector(xgl, 1,ngl);
  free_dvector(wgl, 1,ngl);
  free_imatrix(node, 1,nx, 1,ny);
    
  return 0;  
}//end function init_hex
