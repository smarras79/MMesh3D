/*******************************************************************************
 *
 * SurfMeshParams()
 * 
 * This function comes from Algorithm A9.3 of 
 * [1] L. Piegl and W. Tiller, "The NURBS book", 2nd Ed., 1997. Springe-Verlag
 *
 * It computes the parameters uk and vl to be used to compute the tangent vector
 * at every node Q
 *
 * Input:  n,m, Q[1:q][1:p]
 * Output: uk,vl
 * 
 * Q     = set of data points (those to be interpolated through the surface defined by P).
 * uk,vl = parameters used in the definition of additional quantities during the interp. process.
 *
 * Sept. 11, 2010
 * simone.marras@gmail.com
 *
 *******************************************************************************/

void SURF_MESH_PARAMS(int n, int m, double **Q, double *uk, double *vl)
{
  int num,total;
  int k,l;

  float d;

  /***********************************************
   *  First get the uk:
   ***********************************************/
 
  //Number of nondegenerate rows
  num = m + 1;
  
  //Compute uk[k]:
  uk[0] = 0.0;
  uk[n] = 1.0;

  for(k=1; k<n; k++)
    uk[k] = 0.0;

  for(l=0; l<=m; l++)
  {
    //Total chord length of row
    total = 0.0;
    for(k=1; k<=n; k++)
      {
	cds[k] = Distance3D(Q[k][l][1],Q[K][l][2],Q[k][l][3], Q[k-1][l][1],Q[k-1][l][2],Q[k-1][l][3]);
	total = total + cds[k];
      }
    if(total == 0.0)
      num = num - 1;
    else
      {
	d = 0.0;
	for(k=1; k<n; k++)
	  {
	    d = d + cds[k];
	    uk[k] = uk[k] + d/total;
	  }
      }
  }
  if(num == 0) return(error);
  for(k=1; k<n; k++)
    {
      uk[k] = uk[k]/num;
    }

  //Compute vl[l]:
  vl[0] = 0.0;
  vl[m] = 1.0;

  for(l=1; l<m; l++)
    vl[l] = 0.0;

  for(l=0; l<=m; l++)
  {
    //Total chord length of row
    total = 0.0;
    for(l=1; l<=m; l++)
      {
	cds[l] = Distance3D(Q[k][l][1],Q[K][l][2],Q[k][l][3], Q[k][l-1][1],Q[k][l-1][2],Q[k][l-1][3]);
	total = total + cds[l];
      }
    if(total == 0.0)
      num = num - 1;
    else
      {
	d = 0.0;
	for(l=1; l<m; l++)
	  {
	    d = d + cds[l];
	    vl[l] = vl[l] + d/total;
	  }
      }
  }
  if(num == 0) return(error);
  for(l=1; l<m; l++)
    {
      vl[l] = vl[l]/num;
    }



  return;
}

/*******************************************************************************
 *
 * Distance3D()
 * 
 * This function computes the distance between two 3D points
 *
 * Sept. 11, 2010
 * simone.marras@gmail.com
 *
 *******************************************************************************/

float Distance3D(float pt1x, float pt1y, float pt1z, float pt2x, float pt2y, float pt2z)
{
  float d;

  d = fsqrt( (pt1x - pt2x)*(pt1x - pt2x) + \ //x1-x2
	     (pt1y - pt2y)*(pt1y - pt2y) + \ //y1-y2
	     (pt1z - pt2z)*(pt1z - pt2z));   //z1-z2

  return d;
}


/*******************************************************************************
 *
 * local_surf_interp()
 * 
 * This function comes from Algorithm A9.5 of 
 * [1] L. Piegl and W. Tiller, "The NURBS book", 2nd Ed., 1997. Springe-Verlag
 *
 * Local surface interpolation through (n+1)x(m+1) given points.
 *
 * Input:  n,m, Q[1:q][1:p]
 * Output: U,V,P[1:n][1:m]
 * 
 * Q   = set of data points (those to be interpolated through the surface defined by P).
 * P   = interpolating surface coordinates.
 * n,m = number of points in x and y to be used as interpolated nodes. 
 * In general (n,m) > number of data nodes (q,p) from Q
 *
 * Sept. 11, 2010
 * simone.marras@gmail.com
 *
 *******************************************************************************/

void LOCAL_SURF_INTERP(int n, int m, double **Q, double *U, double *V, double **P)
{
  int total;
  int i,j,k;

  float ub[n+1],vb[m+1];
  float td[n+1][m+1][3];
  float r[m+1],s[n+1];    //Total chord lengths of the rows (r) and columns (s)
  float duk,dvl;
  float alpha, qk;

  //Compute uk, vl:
  SURF_MESH_PARAMS(n, m, Q,  ub, vb);


  //Get ub[], r[] and u direction tangents
  total = 0.0;
  for(k=0; k<=n; k++) 
    ub[k] = 0.0;
  
  for(l=0; l<=m; l++)
    {
      //Compute and load T0,l into td[0][l][0]
      if(l == 0)
	{
	  alpha = 1.0;
	  qk = Q[][][]
	}
 duk = ub[k] - 
    
	r[l] = 0.0;
    }

  return;
}
