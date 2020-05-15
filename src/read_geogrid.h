/* File: read_geogrid.h */

int read_geogrid(
      char * fname,            /* The name of the file to read from */
      int * len,               /* The length of the filename */
      float * rarray,          /* The array to be filled */
      int * nx,                /* x-dimension of the array */
      int * ny,                /* y-dimension of the array */
      int * nz,                /* z-dimension of the array */
      int * isigned,           /* 0=unsigned data, 1=signed data */
      int * endian,            /* 0=big endian, 1=little endian */
      float * scalefactor,     /* value to multiply array elements by before truncation to integers */
      int * wordsize,          /* number of bytes to use for each array element */
      int * status);

void READ_INDEX(char *input_inp, char **topo_type, char **topo_sign, char **topo_projection, char **topo_units, char **topo_description, double *dx, double *dy, double *known_x, double *known_y, double *known_lat, double *known_lon, int *wordsize, int *tile_x, int *tile_y, int *tile_z, int *tile_bdr);
