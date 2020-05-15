/*READ_INPUT.h*/

//int READ_INPUT(char *input_file, int *INPUTVariables, double **BDY_COORDS, double *parameters, char *problem[]);
int READ_INPUT(char *input_file);

int read_parameter_file(double *parameters, char *vertical_coords[]);

int CountNRows(char *input_file);

void load_input_fvect(char *file_name, float *input_vect);

void load_input_ivect(char *file_name, int *input_vect);
