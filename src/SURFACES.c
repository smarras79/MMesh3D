/**********************************************************************
 *
 * SURFACES.c
 *
 * This function builds volume grid:
 *
 **********************************************************************/

#include "myinclude.h"
#include "mydefine.h"

int SURFACES(int isurface, int nnodesu, int nnodesv, int nnodesw, int nelem, \
	     double **BDY_COORDS, char *problem[], double *parameters, double **COORDS, double *topo_surface, double deltaLon, double deltaLat)
{

  int u,v,w,ipoin,i,j,k;
  double *eps,*eta,*xsi;

  unsigned int sigma_flg, hybrid_flg, sleve_flg, tfi_flg;

  double *x,*y,*z;
  double x1x,x2x,x3x,x12x,x13x,x23x,x123x;
  double x1y,x2y,x3y,x12y,x13y,x23y,x123y;
  double x1z,x2z,x3z,x12z,x13z,x23z,x123z;

  double xmin,xmax,ymin,ymax,zmin,zmax;
  double xlength,ylength,zlength;
  double dx,dy,dz;
  double zij,zijminus1;

  //6 boundary surfaces:
  float ***Rbott, ***Rtop;             //Surfaces: R[1:2][1:u][1:v]
  float ***Rwest, ***Reast;            //Surfaces: R[1:2][1:v][1:w]
  float ***Rnorth,***Rsouth;           //Surfaces: R[1:2][1:u][1:w]
  
  //12 boundary edges:
  double **xsw, **xse, **xnw, **xne;
  double **xbw, **xtw, **xbe, **xte;
  double **xbs, **xts, **xbn, **xtn;
  
  //8 corner points:
  double *xwsb, *xwst, *xwnb, *xwnt;
  double *xesb, *xest, *xenb, *xent;

  //Vertical coordinate variables:
  double *sigma;

  //Dummy variables used as dummy mid points in using the function parabola.c
  double midpt_x, midpt_y, midpt_z;   
 
  //Coefficientis of the parabolic distributions.Computed in parabola.c
  double H,ax,ay,az,bx,by,bz,cx,cy,cz;

  double h0,a,b,b1,b2,s1,s2,h1,h2;
  //  double *parameters;
  char *vertical_coords[1];

  //Bdy flag file: instead of sotirng the boundaries in an array, I write the flag 
  //to the file at every loop
  int BDYFLAG[4];
  FILE *BDYFLAGf_ID; 
  char *outfile_bdyflag;
  char nodesu[5], nodesv[5], nodesw[5];
  outfile_bdyflag = (char*) malloc(96 * sizeof(char *));

  //Bdy flag output file:
  sprintf(nodesu, "%d", nnodesu);
  sprintf(nodesv, "%d", nnodesv);
  sprintf(nodesw, "%d", nnodesw);
  
  strcpy(outfile_bdyflag, "3D_mesh_");
  strcat(outfile_bdyflag, problem[2]);
  strcat(outfile_bdyflag,"_");
  strcat(outfile_bdyflag,nodesu);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesv);
  strcat(outfile_bdyflag,"x");
  strcat(outfile_bdyflag,nodesw);
  strcat(outfile_bdyflag, "_BDYFLAG.dat");

  /*****************************************************************************
   * Dynamic allocation where needed:
   *****************************************************************************/
  x = dvector(1,nnodesu);
  y = dvector(1,nnodesv);
  z = dvector(1,nnodesw);

  /*EPS,ETA,XSI*/
  eps = dvector(1,nnodesu);
  eta = dvector(1,nnodesv);
  xsi = dvector(1,nnodesw);

  /*SIGMA*/
  sigma = dvector(1,nnodesw);

  /*Volume surfaces*/
  Rtop   = f3tensor(1,nnodesu, 1,nnodesv, 1,3);
  Rbott  = f3tensor(1,nnodesu, 1,nnodesv, 1,3);
  Reast  = f3tensor(1,nnodesv, 1,nnodesw, 1,3);
  Rwest  = f3tensor(1,nnodesv, 1,nnodesw, 1,3);
  Rsouth = f3tensor(1,nnodesu, 1,nnodesw, 1,3);
  Rnorth = f3tensor(1,nnodesu, 1,nnodesw, 1,3);

  /*Boundary edges*/
  xsw=dmatrix(1,nnodesw,1,3), xse=dmatrix(1,nnodesw,1,3);
  xnw=dmatrix(1,nnodesw,1,3), xne=dmatrix(1,nnodesw,1,3);
  xbw=dmatrix(1,nnodesv,1,3), xtw=dmatrix(1,nnodesv,1,3);
  xbe=dmatrix(1,nnodesv,1,3), xte=dmatrix(1,nnodesv,1,3);
  xbs=dmatrix(1,nnodesu,1,3), xts=dmatrix(1,nnodesu,1,3);
  xbn=dmatrix(1,nnodesu,1,3), xtn=dmatrix(1,nnodesu,1,3);

  /*Volume corner points*/
  xwsb = dvector(1,3), xwst = dvector(1,3), xwnb = dvector(1,3), xwnt = dvector(1,3);
  xesb = dvector(1,3), xest = dvector(1,3), xenb = dvector(1,3), xent = dvector(1,3);

  vertical_coords[0] = (char*) malloc(32 * sizeof(char *));
  //  parameters = dvector(1,NUMBER_OF_PARAMETERS);

  //WARNING: Using this definitions of xmin,etc. implies that the surface that is being
  //         defined has straigth and parallel edges:

  /******************************************************
   *
   * Define the four corners of the domain to be plotted
   *
   * In the case of text file coming from NOAA, we take the
   * simply transform into meters the dimension obtained
   * from the lon lat given in degrees.
   * For every other case (FUNCTION, or NCDF, etc.) the 
   * user input is what counts.
   *
   * Only zmin and zmax are taken from the input in all cases
   *
   ******************************************************/
  if( !strcmp(problem[9], "text")   ||		\
      !strcmp(problem[9], "txt")    ||		\
      !strcmp(problem[9], "noaa")   ||		\
      !strcmp(problem[9], "NOAA")   ||          \
      !strcmp(problem[9], "DEM")   ||		\
      !strcmp(problem[9], "Dem")   ||		\
      !strcmp(problem[9], "dem"))
    {

      //Transform from degrees to meters
      deg2meters(deltaLon, deltaLat, &xlength, &ylength);
      
      xmin = 0.0;
      xmax = xmin + xlength;
      BDY_COORDS[1][1] = xmin;
      BDY_COORDS[2][1] = xmax;

      ymin = 0.0;
      ymax = ymin + ylength;
      BDY_COORDS[1][2] = ymin;
      BDY_COORDS[4][2] = ymax;   

      zmin = 0.0;
      zmax = BDY_COORDS[5][3]; //The domain heigh is still from the input coordinates
      BDY_COORDS[1][3] = zmin;
      BDY_COORDS[2][3] = zmin;
      BDY_COORDS[3][3] = zmin;
      BDY_COORDS[4][3] = zmin;
      BDY_COORDS[6][3] = zmax;
      BDY_COORDS[7][3] = zmax;
      BDY_COORDS[8][3] = zmax;
      zlength = zmax - zmin;

      //Build the coordinates according to the input file
      //and not to the user input coordinates (except for zmax):
      BDY_COORDS[2][2] = BDY_COORDS[1][2];

      BDY_COORDS[3][1] = BDY_COORDS[2][1];
      BDY_COORDS[3][2] = BDY_COORDS[4][2];

      BDY_COORDS[4][1] = xmin;

      BDY_COORDS[5][1] = BDY_COORDS[1][1];
      BDY_COORDS[5][2] = BDY_COORDS[1][2];
    
      BDY_COORDS[6][1] = BDY_COORDS[2][1];
      BDY_COORDS[6][2] = BDY_COORDS[2][2];
     
      BDY_COORDS[7][1] = BDY_COORDS[3][1];
      BDY_COORDS[7][2] = BDY_COORDS[3][2];
     
      BDY_COORDS[8][1] = BDY_COORDS[4][1];
      BDY_COORDS[8][2] = BDY_COORDS[4][2];

      /*printf(" COORDS NOW:\n"); 
      printf(" %f %f %f\n", BDY_COORDS[1][1],BDY_COORDS[1][2],BDY_COORDS[1][3]);
      printf(" %f %f %f\n", BDY_COORDS[2][1],BDY_COORDS[2][2],BDY_COORDS[2][3]);
      printf(" %f %f %f\n", BDY_COORDS[3][1],BDY_COORDS[3][2],BDY_COORDS[3][3]);
      printf(" %f %f %f\n", BDY_COORDS[4][1],BDY_COORDS[4][2],BDY_COORDS[4][3]);
      printf(" %f %f %f\n", BDY_COORDS[5][1],BDY_COORDS[5][2],BDY_COORDS[5][3]);
      printf(" %f %f %f\n", BDY_COORDS[6][1],BDY_COORDS[6][2],BDY_COORDS[6][3]);
      printf(" %f %f %f\n", BDY_COORDS[7][1],BDY_COORDS[7][2],BDY_COORDS[7][3]);
      printf(" %f %f %f\n", BDY_COORDS[8][1],BDY_COORDS[8][2],BDY_COORDS[8][3]);
      */
    }
  else
    {
      xmin = BDY_COORDS[1][1];
      xmax = BDY_COORDS[2][1];
      xlength = xmax - xmin;
      ymin = BDY_COORDS[1][2];
      ymax = BDY_COORDS[4][2];
      ylength = ymax - ymin;
      zmin = BDY_COORDS[1][3];
      zmax = BDY_COORDS[5][3];
      zlength = zmax - zmin;
    }
  
  dlinspace(xmin, xmax, nnodesu, x, "f");
  dlinspace(ymin, ymax, nnodesv, y, "f");
  dlinspace(zmin, zmax, nnodesw, z, "f");

  H = zmax;
  
  /**********************************************************************************
   * Store the parameters of the mountain defined in the TOPOGRAPHY section in Input
   **********************************************************************************/
  //read_parameter_file(parameters, vertical_coords);
                                                    
  h0 = parameters[4];
  a  = parameters[5];
  b  = parameters[9];
  s1 = parameters[7];
  s2 = parameters[8];

  if(s1 == 0.0 )
    s1 = H;
  if(s2 == 0.0 )
    s2 = H;

  //Define the flag for the vertical discretiation coordinate:
  tfi_flg = 1, sigma_flg = 0, sleve_flg = 0, hybrid_flg = 0;
  if( !strcmp(problem[1],"TFI") || !strcmp(problem[1],"tfi") )
    tfi_flg = 1, sigma_flg = 0, sleve_flg = 0, hybrid_flg = 0;
  else if( !strcmp(problem[1],"elliptic") || !strcmp(problem[1],"ELLIPTIC") )
    tfi_flg = 1, sigma_flg = 0, sleve_flg = 0, hybrid_flg = 0;
  else if( !strcmp(problem[1],"SIGMA") || !strcmp(problem[1],"sigma") )
    tfi_flg = 0, sigma_flg = 1, sleve_flg = 0, hybrid_flg = 0;
  else if( !strcmp(problem[1],"SLEVE") || !strcmp(problem[1],"sleve") )
    tfi_flg = 0, sigma_flg = 0, sleve_flg = 1, hybrid_flg = 0;
  else if( !strcmp(problem[1],"HYBRID") || !strcmp(problem[1],"hybrid") )
    tfi_flg = 0, sigma_flg = 0, sleve_flg = 0, hybrid_flg = 1;
  else
    {
      printf(" !!! WARNING Vertical COORDINATE: \n");
      printf("     You did not define a correct vertical coordinate system\n");
      printf("     open 'mountain_parameters' and use one of the following:\n");
      printf("     - SIGMA\n");
      printf("     - HYBRID\n");
      printf("     - SLEVE\n");
      printf("     - TFI\n");
      printf("  By default the TFI will be used\n");
      tfi_flg = 1, sigma_flg = 0, sleve_flg = 0, hybrid_flg = 0;
    }
   
  /*****************************************************************************
   * Define the 8 boundary nodes:
   *****************************************************************************/
  //1] Xwsb
  xwsb[1] = BDY_COORDS[1][1];
  xwsb[2] = BDY_COORDS[1][2];
  xwsb[3] = BDY_COORDS[1][3];
   //2] Xesb
  xesb[1] = BDY_COORDS[2][1];
  xesb[2] = BDY_COORDS[2][2];
  xesb[3] = BDY_COORDS[2][3];
 
  //3] Xenb
  xenb[1] = BDY_COORDS[3][1];
  xenb[2] = BDY_COORDS[3][2];
  xenb[3] = BDY_COORDS[3][3];
  //4] Xwnb
  xwnb[1] = BDY_COORDS[4][1];
  xwnb[2] = BDY_COORDS[4][2];
  xwnb[3] = BDY_COORDS[4][3];
  
  //5] Xwst
  xwst[1] = BDY_COORDS[5][1];
  xwst[2] = BDY_COORDS[5][2];
  xwst[3] = BDY_COORDS[5][3];
  //6] Xest
  xest[1] = BDY_COORDS[6][1];
  xest[2] = BDY_COORDS[6][2];
  xest[3] = BDY_COORDS[6][3];
  
  //7] Xenb
  xent[1] = BDY_COORDS[7][1];
  xent[2] = BDY_COORDS[7][2];
  xent[3] = BDY_COORDS[7][3];
  //8] Xwnb
  xwnt[1] = BDY_COORDS[8][1];
  xwnt[2] = BDY_COORDS[8][2];
  xwnt[3] = BDY_COORDS[8][3];
  
  
  /*****************************************************************************
   * Build the 12 boundary edges of the volume:
   *****************************************************************************/

  //1] Xbs:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwsb, xesb, "y(x)", &ay, &by, &cy);
  parabola2pts(xwsb, xesb, "z(x)", &az, &bz, &cz);

  //Divide the BS side:
  xbs[1][1] = xwsb[1]; //x of Xwsb
  xbs[1][2] = xwsb[2]; //y of Xwsb
  xbs[1][3] = xwsb[3]; //z of Xwsb

  dx = (xesb[1] - xwsb[1])/(nnodesu-1);
  for (u=2; u<=nnodesu; u++)
    {
      xbs[u][1] = xbs[u-1][1] + dx;
      xbs[u][2] = ay*(xbs[u][1])*(xbs[u][1]) + by*xbs[u][1] + cy;
      xbs[u][3] = az*(xbs[u][1])*(xbs[u][1]) + bz*xbs[u][1] + cz;
    }
  
  //2] Xbn:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwnb, xenb, "y(x)", &ay, &by, &cy);
  parabola2pts(xwnb, xenb, "z(x)", &az, &bz, &cz);

  //Divide the BS side:
  xbn[1][1] = xwnb[1]; //x of Xwnb
  xbn[1][2] = xwnb[2]; //y of Xwnb
  xbn[1][3] = xwnb[3]; //z of Xwnb

  dx = (xenb[1] - xwnb[1])/(nnodesu-1);
  for (u=2; u<=nnodesu; u++)
    {
      xbn[u][1] = xbn[u-1][1] + dx;
      xbn[u][2] = ay*(xbn[u][1])*(xbn[u][1]) + by*xbn[u][1] + cy;
      xbn[u][3] = az*(xbn[u][1])*(xbn[u][1]) + bz*xbn[u][1] + cz;
    }

  //3] Xtn:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwnt, xent, "y(x)", &ay, &by, &cy);
  parabola2pts(xwnt, xent, "z(x)", &az, &bz, &cz);

  //Divide the TSN side:
  xtn[1][1] = xwnt[1]; //x of Xwnb
  xtn[1][2] = xwnt[2]; //y of Xwnb
  xtn[1][3] = xwnt[3]; //z of Xwnb

  dx = (xent[1] - xwnt[1])/(nnodesu-1);
  for (u=2; u<=nnodesu; u++)
    {
      xtn[u][1] = xtn[u-1][1] + dx;
      xtn[u][2] = ay*(xtn[u][1])*(xtn[u][1]) + by*xtn[u][1] + cy;
      xtn[u][3] = az*(xtn[u][1])*(xtn[u][1]) + bz*xtn[u][1] + cz;
    }

  //4] Xts:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwst, xest, "y(x)", &ay, &by, &cy);
  parabola2pts(xwst, xest, "z(x)", &az, &bz, &cz);

  //Divide the TS side:
  xts[1][1] = xwst[1]; //x of Xwsb
  xts[1][2] = xwst[2]; //y of Xwsb
  xts[1][3] = xwst[3]; //z of Xwsb

  dx = (xest[1] - xwst[1])/(nnodesu-1);
  for (u=2; u<=nnodesu; u++)
    {
      xts[u][1] = xts[u-1][1] + dx;
      xts[u][2] = ay*(xts[u][1])*(xts[u][1]) + by*xts[u][1] + cy;
      xts[u][3] = az*(xts[u][1])*(xts[u][1]) + bz*xts[u][1] + cz;
    }
 
  /**************************************************************/

  //5] Xbe:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xesb, xenb, "x(y)", &ax, &bx, &cx);
  parabola2pts(xesb, xenb, "z(y)", &az, &bz, &cz);

  //Divide the TE side:
  xbe[1][1] = xesb[1]; //x of Xesb
  xbe[1][2] = xesb[2]; //y of Xesb
  xbe[1][3] = xesb[3]; //z of Xesb

  dy = (xenb[2] - xesb[2])/(nnodesv-1);
  for (v=2; v<=nnodesv; v++)
    {
      xbe[v][2] = xbe[v-1][2] + dy;
      xbe[v][1] = ax*(xbe[v][2])*(xbe[v][2]) + bx*xbe[v][2] + cx;
      xbe[v][3] = az*(xbe[v][2])*(xbe[v][2]) + bz*xbe[v][2] + cz;
    }

  //6] Xte:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xest, xent, "x(y)", &ax, &bx, &cx);
  parabola2pts(xest, xent, "z(y)", &az, &bz, &cz);

  //Divide the BS side:
  xte[1][1] = xest[1]; //x of Xest
  xte[1][2] = xest[2]; //y of Xest
  xte[1][3] = xest[3]; //z of Xest

  dy = (xent[2] - xest[2])/(nnodesv-1);
  for (v=2; v<=nnodesv; v++)
    {
      xte[v][2] = xte[v-1][2] + dy;
      xte[v][1] = ax*(xte[v][2])*(xte[v][2]) + bx*xte[v][2] + cx;
      xte[v][3] = az*(xte[v][2])*(xte[v][2]) + bz*xte[v][2] + cz;
    }

  //7] Xbw:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwsb, xwnb, "x(y)", &ax, &bx, &cx);
  parabola2pts(xwsb, xwnb, "z(y)", &az, &bz, &cz);

  //Divide the BW side:
  xbw[1][1] = xwsb[1]; //x of Xwsb
  xbw[1][2] = xwsb[2]; //y of Xwsb
  xbw[1][3] = xwsb[3]; //z of Xwsb

  dy = (xwnb[2] - xwsb[2])/(nnodesv-1);
  for (v=2; v<=nnodesv; v++)
    {
      xbw[v][2] = xbw[v-1][2] + dy;
      xbw[v][1] = ax*(xbw[v][2])*(xbw[v][2]) + bx*xbw[v][2] + cx;
      xbw[v][3] = az*(xbw[v][2])*(xbw[v][2]) + bz*xbw[v][2] + cz;
    }
  
  //8] XTW:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwst, xwnt, "x(y)", &ax, &bx, &cx);
  parabola2pts(xwst, xwnt, "z(y)", &az, &bz, &cz);

  //Divide the BS side:
  xtw[1][1] = xwst[1]; //x of Xest
  xtw[1][2] = xwst[2]; //y of Xest
  xtw[1][3] = xwst[3]; //z of Xest

  dy = (xwnt[2] - xwst[2])/(nnodesv-1);
  for (v=2; v<=nnodesv; v++)
    {
      xtw[v][2] = xtw[v-1][2] + dy;
      xtw[v][1] = ax*(xtw[v][2])*(xtw[v][2]) + bx*xtw[v][2] + cx;
      xtw[v][3] = az*(xtw[v][2])*(xtw[v][2]) + bz*xtw[v][2] + cz;
    }
  
  /**************************************************************/
  
  //9] XSW:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwsb, xwst, "x(z)", &ax, &bx, &cx);
  parabola2pts(xwsb, xwst, "y(z)", &ay, &by, &cy);

  //Divide the BS side:
  xsw[1][1] = xwsb[1]; //x of Xest
  xsw[1][2] = xwsb[2]; //y of Xest
  xsw[1][3] = xwsb[3]; //z of Xest
  
  w=0;
  dz = (xwst[3] - xwsb[3])/(nnodesw-1);
  for (w=2; w<=nnodesw; w++)
    {
      xsw[w][3] = xsw[w-1][3] + dz;
      xsw[w][1] = ax*(xsw[w][3])*(xsw[w][3]) + bx*xsw[w][3] + cx;
      xsw[w][2] = ay*(xsw[w][3])*(xsw[w][3]) + by*xsw[w][3] + cy;
    }

  //10] XNW:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xwnb, xwnt, "x(z)", &ax, &bx, &cx);
  parabola2pts(xwnb, xwnt, "y(z)", &ay, &by, &cy);

  //Divide the BS side:
  xnw[1][1] = xwnb[1]; //x of Xwnb
  xnw[1][2] = xwnb[2]; //y of Xwnb
  xnw[1][3] = xwnb[3]; //z of Xwnb

  dz = (xwnt[3] - xwnb[3])/(nnodesw-1);
  for (w=2; w<=nnodesw; w++)
    {
      xnw[w][3] = xnw[w-1][3] + dz;
      xnw[w][1] = ax*(xnw[w][3])*(xnw[w][3]) + bx*xnw[w][3] + cx;
      xnw[w][2] = ay*(xnw[w][3])*(xnw[w][3]) + by*xnw[w][3] + cy;
    }
 
  //11] XSE:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xesb, xest, "x(z)", &ax, &bx, &cx);
  parabola2pts(xesb, xest, "y(z)", &ay, &by, &cy);

  //Divide the BS side:
  xse[1][1] = xesb[1]; //x of Xesb
  xse[1][2] = xesb[2]; //y of Xesb
  xse[1][3] = xesb[3]; //z of Xesb

  dz = (xest[3] - xesb[3])/(nnodesw-1);
  for (w=2; w<=nnodesw; w++)
    {
      xse[w][3] = xse[w-1][3] + dz;
      xse[w][1] = ax*(xse[w][3])*(xse[w][3]) + bx*xse[w][3] + cx;
      xse[w][2] = ay*(xse[w][3])*(xse[w][3]) + by*xse[w][3] + cy;
    }

  //12] XNE:
  //Find the parabola coefficients a,b,c for these three points:
  parabola2pts(xenb, xent, "x(z)", &ax, &bx, &cx);
  parabola2pts(xenb, xent, "y(z)", &ay, &by, &cy);

  //Divide the BS side:
  xne[1][1] = xenb[1]; //x of Xwnb
  xne[1][2] = xenb[2]; //y of Xwnb
  xne[1][3] = xenb[3]; //z of Xwnb

  dz = (xent[3] - xenb[3])/(nnodesw-1);
  for (w=2; w<=nnodesw; w++)
    {
      xne[w][3] = xne[w-1][3] + dz;
      xne[w][1] = ax*(xne[w][3])*(xne[w][3]) + bx*xne[w][3] + cx;
      xne[w][2] = ay*(xne[w][3])*(xne[w][3]) + by*xne[w][3] + cy;
    }

  
  /*****************************************************************************
   * Build the 6 boundary surfaces of the volume:
   *****************************************************************************/
  
  //Bottom surface:
  ipoin = 0;
  if( !strcmp(problem[9], "NCDF")   ||		\
      !strcmp(problem[9], "ncdf")   ||		\
      !strcmp(problem[9], "netcdf") ||		\
      !strcmp(problem[9], "NetCDF") ||		\
      !strcmp(problem[9], "Ncdf")   ||		\
      !strcmp(problem[9], "text")   ||		\
      !strcmp(problem[9], "txt")    ||		\
      !strcmp(problem[9], "noaa")   ||		\
      !strcmp(problem[9], "NOAA")   ||          \
      !strcmp(problem[9], "DEM")    ||          \
      !strcmp(problem[9], "Dem")    ||          \
      !strcmp(problem[9], "dem"))
    {
      for(v=1; v<=nnodesv; v++){
	for(u=1; u<=nnodesu; u++){
  
	  Rbott[u][v][1] = x[u];
	  Rbott[u][v][2] = y[v];
	  Rbott[u][v][3] = topo_surface[ipoin+1];
	  ipoin = ipoin + 1;
	}
      }
    }
  else if( !strncmp(problem[9], "FUNCTION", 4) || \
	   !strncmp(problem[9], "function", 4) || \
	   !strncmp(problem[9], "user", 4))
    {
      for(v=1; v<=nnodesv; v++){
	for(u=1; u<=nnodesu; u++){
	  Rbott[u][v][1] = x[u];
	  Rbott[u][v][2] = y[v];
	  Rbott[u][v][3] = topo_user_3Dfunction(problem[5], x[u], y[v], a, b, h0, 3.0/2.0);
	  ipoin = ipoin + 1;
	}
      }
    }
  
  //Free memory taken by the topo_sirface array (rarray). It is not used ever again from this point on:
  //free_dvector(topo_surface, 0,nnodesu*nnodesv);

  u=0;
  v=0;
  //Top surface:
  for(v=1; v<=nnodesv; v++){
    for(u=1; u<=nnodesu; u++){
      Rtop[u][v][1] = x[u];
      Rtop[u][v][2] = y[v];
      Rtop[u][v][3] = zmax;
      }
  }
  v=0;
  w=0;
  //East surface:
  for(w=1; w<=nnodesw; w++){
    for(v=1; v<=nnodesv; v++){
      Reast[v][w][1] = xmax;
      Reast[v][w][2] = y[v];
      Reast[v][w][3] = z[w];
     }
  }
  v=0;
  w=0;
  //West surface:
  for(w=1; w<=nnodesw; w++){
    for(v=1; v<=nnodesv; v++){
      Rwest[v][w][1] = xmin;
      Rwest[v][w][2] = y[v];
      Rwest[v][w][3] = z[w];
    }
  }
  u=0;
  w=0;
  //South surface:
  for(w=1; w<=nnodesw; w++){
    for(u=1; u<=nnodesu; u++){
      Rsouth[u][w][1] = x[u];
      Rsouth[u][w][2] = ymin;
      Rsouth[u][w][3] = z[w];
    }
  }
  u=0;
  w=0;
  //North surface:
  for(w=1; w<=nnodesw; w++){
    for(u=1; u<=nnodesu; u++){
      Rnorth[u][w][1] = x[u];
      Rnorth[u][w][2] = ymax;
      Rnorth[u][w][3] = z[w];
    }
  }


  /*****************************************************************************
   * Define the vertical coordinate systme: SIGMA, HYBRID, SLEVE
   *****************************************************************************/

  for (u=1; u<=nnodesu; u++)
    {
      for (v=1; v<=nnodesv; v++)
	{
	  zijminus1 = Rbott[u][v][3];
	  for (w=2; w<=nnodesw; w++)
	    {
	      dz = (H-Rbott[u][v][3])/(nnodesw-1);
	      zij = zijminus1 + dz;
	      sigma[w] = H*( zij - Rbott[u][v][3])/(H - Rbott[u][v][3]);
	      zijminus1 = zij;
	    }
	}
    }

  /*****************************************************************************
   * Mesh the surfaces:
   *****************************************************************************/
 
  dlinspace(0,1,nnodesu, eps, "force");
  dlinspace(0,1,nnodesv, eta, "force");
  dlinspace(0,1,nnodesw, xsi, "force");

  ipoin = 0;
  u=0;
  v=0;
  w=0;
  
  //Open the bdy flag file for writing:
  if( !strcmp(problem[12], "yes" ))
    BDYFLAGf_ID = fopen(outfile_bdyflag, "a+");
  
  for (w=1; w<=nnodesw; w++){
    for (v=1; v<=nnodesv; v++){
      for (u=1; u<=nnodesu; u++){
	
	/*DA AGGIUSTARE PERCHE le superfici R non sono ancora definite per x,y,z!!! */
	x1x = (1-eps[u])*Rwest[v][w][1]  + eps[u]*Reast[v][w][1];
	x2x = (1-eta[v])*Rsouth[u][w][1] + eta[v]*Rnorth[u][w][1];
	x3x = (1-xsi[w])*Rbott[u][v][1]  + xsi[w]*Rtop[u][v][1];
	
	x1y = (1-eps[u])*Rwest[v][w][2]  + eps[u]*Reast[v][w][2];
	x2y = (1-eta[v])*Rsouth[u][w][2] + eta[v]*Rnorth[u][w][2];
	x3y = (1-xsi[w])*Rbott[u][v][2]  + xsi[w]*Rtop[u][v][2];

	x1z = (1-eps[u])*Rwest[v][w][3]  + eps[u]*Reast[v][w][3];
	x2z = (1-eta[v])*Rsouth[u][w][3] + eta[v]*Rnorth[u][w][3];
	x3z = (1-xsi[w])*Rbott[u][v][3]  + xsi[w]*Rtop[u][v][3];
	
	x12x = (1-eps[u])*(1-eta[v])*xsw[w][1] + (1-eps[u])*eta[v]*xnw[w][1] + \
	  eps[u]*(1-eta[v])*xse[w][1] + eps[u]*eta[v]*xne[w][1];
	x13x = (1-eps[u])*(1-xsi[w])*xbw[v][1] + (1-eps[u])*xsi[w]*xtw[v][1] + \
	  eps[u]*(1-xsi[w])*xbe[v][1] + eps[u]*xsi[w]*xte[v][1];
	x23x = (1-eta[v])*(1-xsi[w])*xbs[u][1] + (1-eta[v])*xsi[w]*xts[u][1] + \
	  eta[v]*(1-xsi[w])*xbn[u][1] + eta[v]*xsi[w]*xtn[u][1];

	x12y = (1-eps[u])*(1-eta[v])*xsw[w][2] + (1-eps[u])*eta[v]*xnw[w][2] + \
	  eps[u]*(1-eta[v])*xse[w][2] + eps[u]*eta[v]*xne[w][2];
	x13y = (1-eps[u])*(1-xsi[w])*xbw[v][2] + (1-eps[u])*xsi[w]*xtw[v][2] + \
	  eps[u]*(1-xsi[w])*xbe[v][2] + eps[u]*xsi[w]*xte[v][2];
	x23y = (1-eta[v])*(1-xsi[w])*xbs[u][2] + (1-eta[v])*xsi[w]*xts[u][2] + 
	  eta[v]*(1-xsi[w])*xbn[u][2] + eta[v]*xsi[w]*xtn[u][2];

	x12z = (1-eps[u])*(1-eta[v])*xsw[w][3] + (1-eps[u])*eta[v]*xnw[w][3] + \
	  eps[u]*(1-eta[v])*xse[w][3] + eps[u]*eta[v]*xne[w][3];
	x13z = (1-eps[u])*(1-xsi[w])*xbw[v][3] + (1-eps[u])*xsi[w]*xtw[v][3] + \
	  eps[u]*(1-xsi[w])*xbe[v][3] + eps[u]*xsi[w]*xte[v][3];
	x23z = (1-eta[v])*(1-xsi[w])*xbs[u][3] + (1-eta[v])*xsi[w]*xts[u][3] + \
	  eta[v]*(1-xsi[w])*xbn[u][3] + eta[v]*xsi[w]*xtn[u][3];


	x123x = (1-eps[u])*(1-eta[v])*(1-xsi[w])*xwsb[1] + \
	  (1-eps[u])*(1-eta[v])*xsi[w]*xwst[1]           + \
	  (1-eps[u])*eta[v]*(1-xsi[w])*xwnb[1]           + \
	  (1-eps[u])*eta[v]*xsi[w]*xwnt[1]               + \
	  eps[u]*(1-eta[v])*(1-xsi[w])*xesb[1]           + \
	  eps[u]*(1-eta[v])*xsi[w]*xest[1]               + \
	  eps[u]*eta[v]*(1-xsi[w])*xenb[1]               + \
	  eps[u]*eta[v]*xsi[w]*xent[1];

	x123y = (1-eps[u])*(1-eta[v])*(1-xsi[w])*xwsb[2] + \
	  (1-eps[u])*(1-eta[v])*xsi[w]*xwst[2]           + \
	  (1-eps[u])*eta[v]*(1-xsi[w])*xwnb[2]           + \
	  (1-eps[u])*eta[v]*xsi[w]*xwnt[2]               + \
	  eps[u]*(1-eta[v])*(1-xsi[w])*xesb[2]           + \
	  eps[u]*(1-eta[v])*xsi[w]*xest[2]               + \
	  eps[u]*eta[v]*(1-xsi[w])*xenb[2]               + \
	  eps[u]*eta[v]*xsi[w]*xent[2];

	x123z = (1-eps[u])*(1-eta[v])*(1-xsi[w])*xwsb[3] + \
	  (1-eps[u])*(1-eta[v])*xsi[w]*xwst[3]           + \
	  (1-eps[u])*eta[v]*(1-xsi[w])*xwnb[3]           + \
	  (1-eps[u])*eta[v]*xsi[w]*xwnt[3]               + \
	  eps[u]*(1-eta[v])*(1-xsi[w])*xesb[3]           + \
	  eps[u]*(1-eta[v])*xsi[w]*xest[3]               + \
	  eps[u]*eta[v]*(1-xsi[w])*xenb[3]               + \
	  eps[u]*eta[v]*xsi[w]*xent[3];

	ipoin = ipoin + 1;

	if( !strcmp(problem[1],"SLEVE") || !strcmp(problem[1],"sleve") )
	  {
	    if( Rbott[u][v][1] >= -a && Rbott[u][v][1] <= a)
	      h1 = 0.5*h0*cos(PI*Rbott[u][v][1]/(2.0*sqrt(a*b)));
	    else
	      h1 = 0.0;
	
	    h2 = Rbott[u][v][3] - h1;
	    
	    b1 = sinh((H - sigma[w])/s1)/sinh(H/s1);
	    b2 = sinh((H - sigma[w])/s2)/sinh(H/s2);
	  }
	
	COORDS[ipoin][1] = ipoin;
	COORDS[ipoin][2] = x1x + x2x + x3x - x12x - x13x - x23x + x123x;
	COORDS[ipoin][3] = x1y + x2y + x3y - x12y - x13y - x23y + x123y;
	COORDS[ipoin][4] = tfi_flg*(x1z + x2z + x3z - x12z - x13z - x23z + x123z)       \
	  + sigma_flg  *(Rbott[u][v][3] + sigma[w]*(H - Rbott[u][v][3])/H)              \
	  + hybrid_flg *(sigma[w] + Rbott[u][v][3]*sinh((H - sigma[w])/s1)/sinh(H/s1))  \
	  + sleve_flg  *(sigma[w] + h1*b1 + h2*b2);

	//Assign the boundary flags:
	//a) BOTTOM CORNERS:
	if( u == 1 && v == 1 && w == 1){
	  BDYFLAG[0] = 1;				// == 1 means that is it a boundary node on the front-bott-left corner.
	}
	else if( u == nnodesu && v == 1 && w == 1){
	  BDYFLAG[0] = 2;				// == 2 means that is it a boundary node on the front-bott-right corner.
	}
	else if( u == nnodesu && v == nnodesv && w == 1){
	  BDYFLAG[0] = 3;				// == 3 means that is it a boundary node on the back-bott-right corner.
	}
	else if( u == 1 && v == nnodesv && w == 1){
	  BDYFLAG[0] = 4;				// == 4 means that is it a boundary node on the back-bott-left corner.
	}
	//b) TOP CORNERS:
	else if( u == 1 && v == 1 && w == nnodesw){
	  BDYFLAG[0] = 8;				// == 8 means that is it a boundary node on the front-top-left corner.
	}
	else if( u == nnodesu && v == 1 && w == nnodesw){
	  BDYFLAG[0] = 5;				// == 5 means that is it a boundary node on the front-top-right corner.
	}
	else if( u == nnodesu && v == nnodesv && w == nnodesw){
	  BDYFLAG[0] = 6;				// == 6 means that is it a boundary node on the back-top-right corner.
	}
	else if( u == 1 && v == nnodesv && w == nnodesw){
	  BDYFLAG[0] = 7;				// == 7 means that is it a boundary node on the back-top-left corner.
	}

	//c) BOTT-EDGES (horizontal):
	else if( v == 1 && w == 1 && (u>1 && u<nnodesu)){
	  BDYFLAG[0] = 9;				// == 9 means that is it a boundary node on the bott-front edge.
	}
	else if( u == nnodesu && w == 1 && (v>1 && v<nnodesv)){
	  BDYFLAG[0] = 10;				// == 10 means that is it a boundary node on the bottom-right edge.
	}
	else if( v == nnodesv && w == 1 && (u>1 && u<nnodesu)){
	  BDYFLAG[0] = 11;				// == 11 means that is it a boundary node on the bottom-back edge.
	}
	else if( u == 1 && w == 1 && (v>1 && v<nnodesv)){
	  BDYFLAG[0] = 12;				// == 12 means that is it a boundary node on the bottom-left edge.
	}

	//d) TOP-EDGES (horizontal):
	else if( v == 1 && w == nnodesw && (u>1 && u<nnodesu)){
	  BDYFLAG[0] = 17;				// == 17 means that is it a boundary node on the top-front edge.
	}
	else if( u == nnodesu && w == nnodesw && (v>1 && v<nnodesv)){
	  BDYFLAG[0] = 18;				// == 18 means that is it a boundary node on the top-right edge.
	}
	else if( v == nnodesv && w == nnodesw && (u>1 && u<nnodesu)){
	  BDYFLAG[0] = 19;				// == 19 means that is it a boundary node on the top-back edge.
	}
	else if( u == 1 && w == nnodesw && (v>1 && v<nnodesv)){
	  BDYFLAG[0] = 20;				// == 20 means that is it a boundary node on the top-left edge.
	}

	//e) VERTICAL-EDGES:
	else if( u == 1 && v == 1 && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 16;				// == 16 means that is it a boundary node on the front-left vedge.
	}
	else if( u == nnodesu && v == 1 && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 13;				// == 13 means that is it a boundary node on the front-right vedge.
	}
	else if( u == nnodesu && v == nnodesv && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 14;				// == 14 means that is it a boundary node on the back-right vedge.
	}
	else if( u == 1 && v == nnodesv && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 15;				// == 15 means that is it a boundary node on the back-left vedge.
	}

	//f) FACES:
	else if( w == 1 && (u<nnodesu && u>1) && (v<nnodesv && v>1)){
	  BDYFLAG[0] = 21;				// == 21 means that is it a boundary node on the bottom face
	}
	else if( w == nnodesw && (u<nnodesu && u>1) && (v<nnodesv && v>1)){
	  BDYFLAG[0] = 23;				// == 23 means that is it a boundary node on the top face
	}
	else if( u == nnodesu && (v>1 && v<nnodesv) && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 22;				// == 22 means that is it a boundary node on the east face
	}
	else if( u == 1 && (v>1 && v<nnodesv) && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 24;				// == 24 means that is it a boundary node on the west face
	}
	else if( v == 1 && (u>1 && u<nnodesu) && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 25;				// == 25 means that is it a boundary node on the front face
	}
	else if( v == nnodesv && (u>1 && u<nnodesu) && (w>1 && w<nnodesw)){
	  BDYFLAG[0] = 26;				// == 26 means that is it a boundary node on the back face
	}
	else
	  {
	    BDYFLAG[0] = 0;		                // == 0 means that is it NOT a boundary node
	  }
	//End assign bdy flags.
	
	//Print boundary file:
	if( !strcmp(problem[12],"yes"))
	  fprintf(BDYFLAGf_ID, "%d %lf %lf %lf\n", BDYFLAG[0], COORDS[ipoin][2],COORDS[ipoin][3],COORDS[ipoin][4]);

      }
    }
  }
      
  //Close the bdy file after use:
  if( !strcmp(problem[12], "yes" ))
    fclose(BDYFLAGf_ID);
  
  /*****************************************************************************
   * Free memory where needed:
   *****************************************************************************/
  free_dvector(x,1,nnodesu);
  free_dvector(y,1,nnodesv);
  free_dvector(z,1,nnodesw);
  
  free_dvector(eps, 1,nnodesu);
  free_dvector(eta, 1,nnodesv);
  free_dvector(xsi, 1,nnodesw);
  
  free_dvector(sigma, 1,nnodesw);


  /*Volume surfaces*/
  free_f3tensor(Rtop,   1,nnodesu, 1,nnodesv, 1,3);
  free_f3tensor(Rbott,  1,nnodesu, 1,nnodesv, 1,3);
  free_f3tensor(Reast,  1,nnodesv, 1,nnodesw, 1,3);
  free_f3tensor(Rwest,  1,nnodesv, 1,nnodesw, 1,3);
  free_f3tensor(Rsouth, 1,nnodesu, 1,nnodesw, 1,3);
  free_f3tensor(Rnorth, 1,nnodesu, 1,nnodesw, 1,3);
  
  /*Boundary edges*/
  free_dmatrix(xsw,1,nnodesw,1,3), free_dmatrix(xse,1,nnodesw,1,3), free_dmatrix(xnw,1,nnodesw,1,3), free_dmatrix(xne,1,nnodesw,1,3);
  free_dmatrix(xbw,1,nnodesv,1,3), free_dmatrix(xtw,1,nnodesv,1,3), free_dmatrix(xbe,1,nnodesv,1,3), free_dmatrix(xte,1,nnodesv,1,3);
  free_dmatrix(xbs,1,nnodesu,1,3), free_dmatrix(xts,1,nnodesu,1,3), free_dmatrix(xbn,1,nnodesu,1,3), free_dmatrix(xtn,1,nnodesu,1,3);

  /*volume corner points*/
  free_dvector(xwsb,1,3), free_dvector(xwst,1,3), free_dvector(xwnb,1,3), free_dvector(xwnt,1,3);
  free_dvector(xesb,1,3), free_dvector(xest,1,3), free_dvector(xenb,1,3), free_dvector(xent,1,3);

  //  free_dvector(parameters, 1,NUMBER_OF_PARAMETERS);
  free(vertical_coords[0]);
  free(outfile_bdyflag);

  return 0;
}


/*
 * HIGH-ORDER (nop > 1)
 *

int SURFACES_HIGH_ORDER(int isurface, int nnodesu, int nnodesv, int nnodesw, int nop, int ngl, \
			int nelem, double **BDY_COORDS, char *problem[], double *parameters,   \
			double **COORDS, double *topo_surface, double deltaLon, double deltaLat)
{

  int u,v,w,ipoin,i,j,k;
  double *eps,*eta,*xsi;

  unsigned int sigma_flg, hybrid_flg, sleve_flg, tfi_flg;

  double x,y,z;

  int ***inode;
  int *ele_col;

  double xmin,xmax,ymin,ymax,zmin,zmax;
  double xlength,ylength,zlength;
  double dx,dy,dz;
  double zij,zijminus1;

  //Local variables for High-order
  int ie, ies, ii, jj, kk, ip, i1,j1,k1, ih;
  int nglm1, npoin_g, nelem_g, nboun_g;
  int nx, ny, nz;
  int nelx, nely, nelz;
  int ielex, ieley, ielez, nelex, neley, nelez;
  int ilglx, ilgly, ilglz;
  double xstart, ystart, zstart;
  double Lx, Ly, Lz;

  //LGL points and weigth (weight not really needed here, but just in case):
  double *xgl, *wgl;

  //Bdy flag file: instead of sotirng the boundaries in an array, I write the flag 
  //to the file at every loop
  int BDYFLAG[4];
  FILE *BDYFLAGf_ID;

  char *outfile_bdyflag;
  char *vertical_coords[1];
  char nodesu[5], nodesv[5], nodesw[5];

  outfile_bdyflag = (char*) malloc(96 * sizeof(char *));
  vertical_coords[0] = (char*) malloc(32 * sizeof(char *));
  
  //Number of elements:
  nelex = nnodesu-1;
  neley = nnodesv-1;
  nelez = nnodesw-1;

  //Number of nodes for high-order elements (nop>1):
  nglm1 = ngl - 1;
  nnodesu=nelex*nop + 1;
  nnodesv=neley*nop + 1;
  nnodesw=nelez*nop + 1;
  
  //Bdy flag output file:
  sprintf(nodesu, "%d", nnodesu);
  sprintf(nodesv, "%d", nnodesv);
  sprintf(nodesw, "%d", nnodesw);
  
  /*****************************************************************************
   * Dynamic allocation where needed:
   *****************************************************************************
  inode = i3tensor(1,nnodesu,1,nnodesv,1,nnodesw);
  ele_col = ivector(1,nelem);

  /*EPS,ETA,XSI*
  eps = dvector(1,nnodesu);
  eta = dvector(1,nnodesv);
  xsi = dvector(1,nnodesw);

  /*xgl and wgl: LGL points and weights: *
  xgl = dvector(1,ngl);
  wgl = dvector(1,ngl);
  
  /*****************************************************************************
   * Begin insertion of high-order nodes:
   *****************************************************************************
  //Compute the LGL points:
  legendre_gauss_lobatto(ngl, xgl, wgl);

  xmax = 1;
  xmin = 0;
  ymax = 1;
  ymin = 0;
  zmax = 10;
  zmin = 0;

  Lx = (xmax - xmin)/nelex;
  Ly = (ymax - ymin)/neley;
  Lz = (zmax - zmin)/nelez;

  // Generate Coordinates Here 
  ip = 1;
  ii = 0;
  jj = 0;
  kk = 0;

  for( ielez = 1; ielez <= nelez; ielez++)
    {
      zstart = zmin + (ielez - 1)*Lz;
      if (ielez == 1)
  	k1 = 1;
      else
  	k1 = 2;
      
      for( ilglz = k1; ilglz<=ngl; ilglz++)
  	{
  	  kk = kk + 1;
  	  z = zstart + 0.5*Lz*(1 + xgl[ilglz]);
  	  jj = 0;
  	  for( ieley = 1; ieley <= neley; ieley++)
  	    {
  	      ystart = ymin + (ieley - 1)*Ly;
  	      if (ieley == 1)
  		j1 = 1;
  	      else
  		j1 = 2;
	      
  	      for( ilgly = j1; ilgly <= ngl; ilgly++)
  		{
  		  jj = jj + 1;
              
  		  y= ystart + 0.5*Ly*(1 + xgl[ilgly]);
  		  ii = 0;
		  
  		  for( ielex = 1; ielex <= nelex; ielex++)
  		    {
  		      xstart = xmin + (ielex - 1)*Lx;
		      
  		      if (ielex == 1)
  			i1 = 1;
  		      else
  			i1 = 2;
		      
  		      for( ilglx = i1; ilglx <= ngl; ilglx++)
  			{
  			  ii = ii + 1;
  			  x = xstart + 0.5*Lx*(1 + xgl[ilglx]);
  			  ie = ielex + nelex*(ieley -1) + nelex*neley*(ielez -1);
			  COORDS[ip][1] = ip;
  			  COORDS[ip][2] = x;
  			  COORDS[ip][3] = y;
  			  COORDS[ip][4] = z;
  			  inode[ii][jj][kk] = ip;
			  //printf("X,Y,Z: %d, %lf %lf %lf\n", ip, COORDS[ip][2],COORDS[ip][3],COORDS[ip][4]);
  			  ip = ip + 1;
			}
  		    }
  		}
  	    }
  	}
    }
  
  //Construct Intma
    ie = 0;
    for( ielez = 1; ielez<=nelez; ielez++)
      {
	ies = 0;
	for( ieley = 1; ieley<=neley; ieley++)
	  {
	    for( ielex = 1; ielex<=nelex; ielex++)
	      {
		ies = ies + 1;
		ie = ie + 1;
		ele_col[ie] = ies;
		
		for( ilglz = 1; ilglz<=ngl; ilglz++)
		  {
		    kk = (ielez - 1)*nglm1 + ilglz;
		    for( ilgly = 1; ilgly<=ngl; ilgly++)
		      {
			jj = (ieley - 1)*nglm1 + ilgly;
			for( ilglx = 1; ilglx<=ngl; ilglx++)
			  {
			    ii = (ielex - 1)*nglm1 + ilglx;
			    ih=nx*(ieley-1)*nglm1 + nx*(ilgly-1) + ii;
			    ip=inode[ii][jj][kk];
			    
			    //intma[ilglx][ilgly][ilglz][ie]=ip;
			    
			    //printf(" IE %d %d %d %d %d\n", ilglx, ilgly, ilglz, ie, ip);
			  }
		      }
		  }
	      }
	  }
      }
   
    /*****************************************************************************
     * Free memory where needed:
     *****************************************************************************
    free_i3tensor(inode, 1,nnodesu,1,nnodesv,1,nnodesw);
    free_ivector(ele_col, 1,nelem);

    free_dvector(eps, 1,nnodesu);
    free_dvector(eta, 1,nnodesv);
    free_dvector(xsi, 1,nnodesw);
    
    free_dvector(xgl, 1,ngl);
    free_dvector(wgl, 1,ngl);
    
    //  free_dvector(parameters, 1,NUMBER_OF_PARAMETERS);
    free(vertical_coords[0]);
    free(outfile_bdyflag);
    
    return 0;
} /* END high-order */
