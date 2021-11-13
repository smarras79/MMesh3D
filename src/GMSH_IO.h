char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void gmsh_data_read ( char *gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] );
int *gmsh_mesh2d_element_data_example ( int element_num, int element_order );
void gmsh_mesh2d_element_size_example ( int *element_num, int *element_order );
double *gmsh_mesh2d_node_data_example ( int node_num, int node_dim );
void gmsh_mesh2d_node_size_example ( int *node_num, int *node_dim );
void gmsh_mesh1d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh2d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh3d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_size_read ( char *gmsh_filename, int *node_num, int *node_dim,
  int *element_num, int *element_order );
int *i4mat_copy_new ( int m, int n, int a1[] );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
int s_begin ( char *s1, char *s2 );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
double s_to_r8 ( char *s, int *lchar, int *error );
void timestamp ( );
