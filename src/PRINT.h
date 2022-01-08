/*PRINT.h*/

void VIEW_f2DMAT(char *mat_name, float **A, int min_rows, int max_rows, int min_cols, int max_cols);
void VIEW_d2DMAT(char *mat_name, double **A, int min_rows, int max_rows, int min_cols, int max_cols);
void VIEW_i2DMAT(char *mat_name, int **A, int min_rows, int max_rows, int min_cols, int max_cols);

void VIEW_fVECT(char *vect_name, float *V, int min_rows, int max_rows);
void VIEW_dVECT(char *vect_name, double *V, int min_rows, int max_rows);
void VIEW_iVECT(char *vect_name, int *V, int min_rows, int max_rows);

void VIEW_i3DMAT(char *mat_name, int ***A, int min_rows, int max_rows, int min_cols, int max_cols, int min_depth, int max_depth);

void MAT2F_f2DMAT(char *mat_name, float **A, int min_rows, int max_rows, int min_cols, int max_cols);
void MAT2F_d2DMAT(char *mat_name, double **A, int min_rows, int max_rows, int min_cols, int max_cols);
void MAT2F_i2DMAT(char *mat_name, int **A, int min_rows, int max_rows, int min_cols, int max_cols);

void VECT2F_fVECT(char *vect_name, float *V, int min_rows, int max_rows);
void VECT2F_dVECT(char *vect_name, double *V, int min_rows, int max_rows);
void VECT2F_iVECT(char *vect_name, int *V, int min_rows, int max_rows);

void PRINT(float **u, float **v, float **P, int imax, int jmax);

void PRINT_VECT_VTK(int imax, int jmax, int kmax, char *vector_field_name, float **U, float **V, float **W, float *x, float *y, float *z);

void PRINT_GRID_VTK(char *grid_name, int imax, int jmax, int kmax, float *x, float *y, float *z);

void PRINT_PRESS_VTK(int imax, int jmax, char *press_field_name, float **P);

void PRINT_STRUCT_GRID_VTK(char *grid_name, int imax, int jmax, int kmax, float *x, float *y, float *z);

void dPRINT_UNSTRUCT_GRID_VTK(char *grid_name, int array_numb, int coords_ncolumns);

void dGRID2ALYAFILE(char *outfile_msh_alya, double **COORDS, int **CONN, int nnodes, int nelem, int nperiodic_nodes);

void dPRINT_UNSTRUCT_GRID_GMSH(char *grid_name, int array_numb, int coords_ncolumns);

void dPRINT_UNSTRUCT_GRID_CONN(char *grid_name, int array_numb, int coords_ncolumns);

void dwrt2VTK(char *file_name);

void dwrt2GMSH(char *file_name);

void dwrt2CONN(char *file_name);

int PRINT_WELCOME_MESSAGE(void);

int PRINT_INFO(void);

void PRINT_ERROR(char *);
