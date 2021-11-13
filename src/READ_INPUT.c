/*****************************************************************************
Function to read in the user input file for the Cartesian Mesh Generator:
   
   simone.marras@bsc.es
   April 26, 2008

*****************************************************************************/

#include <errno.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "MEMORY.h"
#include "MYDEFINE.h"
#include "GLOBAL_VARS.h"

#define SEPARATOR " \t"

int READ_INPUT(char *input_file)
{

  /****************************************************************************************
   * Assign the variables and file names read in input, to their correspondings in main.c:
   *
   *
   ****************************************************************************************
   *How the input is stored:
   *
   * Numeric variables:
   *  INPUTVariables[0]    //Number of blocks
   *  INPUTVariables[1]    //Number of elements x block 1
   *  INPUTVariables[2]    //Number of elements x block 2
   *  INPUTVariables[3]    //Number of elements y block 1
   *  INPUTVariables[4]    //Number of elements y block 2
   *  INPUTVariables[5]    //Number of elements z block 1
   *  INPUTVariables[6]    //Number of elements z block 2
   *  INPUTVariables[8]    //Order of the elements (nop = 1 means 1st order element)
   *
   *  Up to now INPUTVariables is allocated in main.c (called INPUTVars) with 
   *  MAX_INPUTS = 20. Increase its value in include_.h if needed.
   *
   * Boundary nodes:
   *  BDY_COORDS[1:nbdy_nodes][1:3] -> x, y, z of every boundary node (by default it's 8 nodes for the cube)
   *
   *		
   * Definition strings:
   *   problem[0]           // Problem name: agnesi, schar, airfoil, etc.
   *   problem[1]           // Meshing scheme: mountain, TFI, etc. followed by two parameters s1 and s2 (double) stored in parameters[7:8]
   *    parameters[7,8]     // --- scheme parameters s1, s2 (used for sleve and hybrid) THESE can be left blank
   * 
   *
   *   problem[2]           // Elements type: QUAD, TRI
   *   problem[3]           // periodicity in X: On/Off
   *   problem[4]           // periodicity in Y: On/Off
   *   problem[8]           // periodicity in Z: On/Off
   *   problem[5]           // Topography file (ex. spain4km.nc)
   *   problem[9]           // Topography file TYPE (ex. NCDF, or netcdf, or noaa, or txt, or ...)
   *   problem[10]          // Mountain parameters file needed when the mountain is user-defined (FUNCTION)
   *   problem[11]          // Flag for TOPOGRAPHY, BATHIMETRY, or BOTH
   *   problem[6,7,12,13]         // Output file flags:
   *                           ALYA:     yes/no
   *                           VTK:      yes/no
   *                           GMSH:     yes/no
   *                           BDY_FILE: yes/no
   *
   *	       	!!! Up to now *problem is allocated in main.c (called *problem[]) with 
   *   		    PROB_ENTRIES = 14 Increase its value in myinclude.h if needed.
   *
   *****************************************************************************************/
  
  FILE *file_ID;
  
  unsigned int count;
  int i,j,k=0, file_nlines, line_cntr=0, nbdy_nodes;
  char header[64],header1[32],header2[32];
  char *ptr_str;

  /* Give default values before reading them into from the input file*/
  lMETIS = 0; /* Use METIS for domain decomposition    */
  lCART  = 1; /* Use MPI_CART for domain decomposition */
  
  
  if((file_ID = fopen(input_file, "r")) == NULL){
    printf(" ERROR The file %s could not be open. It may not be in the working directory.\n", input_file);
    printf(" The program will exit now\n");
    exit(1);
  }
  else{
    fscanf(file_ID, "%s\n", header);
    
    //PROBLEM:
    fscanf(file_ID, "%s\n", header);
    fscanf(file_ID, "%s %s\n", header, header); //Problem name
    strcpy(problem[0], header);
    fscanf(file_ID, "%s\n", header); 
    //END_PROBLEM
    
    //PROBLEM_DEFINITION:
    fscanf(file_ID, "%s\n", header);
    count = 0;
    //fscanf(file_ID, "%s %d\n", header, &INPUTVariables[7]); //number of space dimension (nsd)
    //count++;
    fscanf(file_ID, "%s %s %lf %lf\n", header, header, &parameters[7], &parameters[8]);
    strcpy(problem[1], header);                             //Meshing_scheme and scheme parameters s1, s2 (used for sleve and hybrid)
    count++;

    fscanf(file_ID, "%s %s\n", header, header);             //Element types
    strcpy(problem[2], header);
 
    fscanf(file_ID, "%s %d\n", header, &INPUTVariables[0]); //number of blocks
    count++;

    if(INPUTVariables[0] == 1)
      {
	fscanf(file_ID, "%s %d\n", header, &INPUTVariables[1]); //number of elements in x, block 1
	count++;
	fscanf(file_ID, "%s %d\n", header, &INPUTVariables[3]); //number of elements in y, block 1
	count++;
	fscanf(file_ID, "%s %d\n", header, &INPUTVariables[5]); //number of elements in z, block 1
	count++;
	fscanf(file_ID, "%s %d\n", header, &INPUTVariables[8]); //Order of the elements
	count++;
      }
    else if(INPUTVariables[0] > 1)
      {
	printf("\n");
	printf(" ERROR IN INPUT: \n");
	printf("       Currently, the maximum allowed number of blocks is 1. Ask me how to do it.\n");
	printf("       simone.marras@gmail.com \n");
	printf(" The program will exit now. Change your input and re-run the code.\n");
	printf(" \n");
	exit(1);
      }
    
    fscanf(file_ID, "%s %s\n", header, header);
    strcpy(problem[3], header);                             //Periodicity in X: on/off
    count++;
    fscanf(file_ID, "%s %s\n", header, header);
    strcpy(problem[4], header);                             //Periodicity in Y: on/off
    count++;
    fscanf(file_ID, "%s %s\n", header, header);
    strcpy(problem[8], header);                             //Periodicity in Z: on/off
    count++;
    fscanf(file_ID, "%s\n", header);
    //END_PROBLEM_DEFINITION
    
    //BDY_NODES:
    fscanf(file_ID, "%s\n", header);
    fscanf(file_ID, "%s %s %s %s\n", header, header, header, header);

    //Start storing the boundary node coordinates x,y,z
    MEMORY_ALLOCATE(2); //allocate BDY_COORDS
    
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[1][1], &BDY_COORDS[1][2], &BDY_COORDS[1][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[2][1], &BDY_COORDS[2][2], &BDY_COORDS[2][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[3][1], &BDY_COORDS[3][2], &BDY_COORDS[3][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[4][1], &BDY_COORDS[4][2], &BDY_COORDS[4][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[5][1], &BDY_COORDS[5][2], &BDY_COORDS[5][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[6][1], &BDY_COORDS[6][2], &BDY_COORDS[6][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[7][1], &BDY_COORDS[7][2], &BDY_COORDS[7][3]);
    fscanf(file_ID, "%s %lf %lf %lf\n", header, &BDY_COORDS[8][1], &BDY_COORDS[8][2], &BDY_COORDS[8][3]);
    fscanf(file_ID, "%s\n", header);
    count++;
    //END_BDY_NODES
    
    //TOPOGRAPHY:
    fscanf(file_ID, "%s\n", header);
    fscanf(file_ID, "%s %s %s %s\n", header, header1, header2, header);
    strcpy(problem[9], header1);                           //Topography: TYPE (ex. NCDF, netcdf, FUNCTION/USER, )  
    strcpy(problem[5], header2);                           //If the first word denotes the format of a file, then a file name follows
                                                           //otherwise, if the first word is FUNCTION then a number denoting the function follows.
                                                           //Ex.1 (using a ncdf file to read the topography from):
                                                           //   TOPOGRAPHY:   netcdf  spain4km.nc
                                                           //Ex.2 (using a user-defined function defined in topo_user_function.c
                                                           //   TOPOGRAPHY:   function  1   (1 indicates what case is used in topo_user_function.c)
    strcpy(problem[11], header);                           //TOPOGRAPHY, BATHIMETRY, BOTH, NONE (Flag to define whether you are reading 
                                                           //the negative or positive heigths in the noaa topography file.
    fscanf(file_ID, "%s %lf\n", header, &parameters[4]);   //hm or r
    fscanf(file_ID, "%s %lf\n", header, &parameters[5]);   //a_c, or xc
    fscanf(file_ID, "%s %lf\n", header, &parameters[9]);   //b_c, or yc
    fscanf(file_ID, "%s %lf\n", header, &parameters[6]);   //lambda or yc
    fscanf(file_ID, "%s\n", header);
    count++;
    //END_TOPOGRAPHY

    //OUTPUT_FILES:
    fscanf(file_ID, "%s\n", header);
    fscanf(file_ID, "%s %s\n", header, header); //ALYA
    strcpy(problem[6], header);
    fscanf(file_ID, "%s %s\n", header, header); //VTK
    strcpy(problem[7], header);
    fscanf(file_ID, "%s %s\n", header, header); //GMSH
    strcpy(problem[13], header);
    fscanf(file_ID, "%s %s\n", header, header); //BDY_FILE
    strcpy(problem[12], header);
    fscanf(file_ID, "%s\n", header); 
    //END_PROBLEM
    
    //END_INPUT_FILE
    fscanf(file_ID, "%s\n", header);
  }//End else.
  
  //Close file	
  fclose(file_ID);
  
  // NOTICE on allowed number of variables:
  // Problem[] is now allocated for PROB_ENTRIES 25 defined in include_.h
  // INPUTVarsiables is allocated for MAX_INPUT_ENTRIES 35 idefined n include_.h
  // *varnames[] may contain 6 strings (see its allocation in main.c)
  // If you add more entires to any of these three arrays you need to change the
  // allocation space in include_.h for "*problem[]" and "*varnames[]", and in main.c for "INPUTVariables[]"
  
  /*************************************************************************************
   * Store the input values to global variables:
   *************************************************************************************/
  nelx = INPUTVariables[1];
  nely = INPUTVariables[3];
  nelz = INPUTVariables[5];
  nop  = INPUTVariables[8];
  
  /* For now, the polynomial order is set to be the same in all directions */
  nopx = nop;
  nopy = nop;
  nopz = nop;
  
  ngl  = nop + 1;
  nglx = ngl;
  ngly = ngl;
  nglz = ngl;
  
  nnodesx = nelx*nopx + 1;
  nnodesy = nely*nopy + 1;
  nnodesz = nelz*nopz + 1;
  
  if( !strncmp( problem[2], "HEXA", 4) ){
    EL_NODES = HEXA;
    nelem = nelx*nely*nelz;
  }
  if( !strncmp( problem[2], "WEDGE", 4) ){
    EL_NODES = WEDGE;
    nelem = 2*nelx*nely*nelz;
  }
  else{
    EL_NODES = HEXA;
    nelem = nelx*nely*nelz;
  }
  
  if(nnodesz < 2 ){
    printf(" #------------------------------------------------------------------#\n");
    printf("   !!! ERROR in %s: nnodesz must be at least = 2 !!! \n", inputfile);
    printf("   !!! open %s and set nnodesz = 2 or larger\n", inputfile);
    printf("   !!! The program will exit now\n");
    printf(" #------------------------------------------------------------------#\n");
    exit(1);
  }

  strcpy(print_alya, problem[6]);
  strcpy(print_vtk,  problem[7]);
  strcpy(print_gmsh, problem[13]);
  
  xmin = BDY_COORDS[1][1];
  xmax = BDY_COORDS[2][1];
  ymin = BDY_COORDS[1][2];
  ymax = BDY_COORDS[3][2];
  zmin = BDY_COORDS[1][3];
  zmax = BDY_COORDS[8][3];
  
  xc = (xmin + xmax)/2.0;
  yc = (ymin + ymax)/2.0;
  
  return 0;
}; //END Read input file.


/*******************************************************************
 * read_parameter_file.c
 *
 * this function reads a text file (exisitng within the working directory)
 * with the decay parameters used in these z-coordinate mountain functions
 *
 * Used by: 
 *			mountain_transformations.c
 *			external_geometries.c
 *
 * The file has this structure:
 *
 * hm:				450.0
 * a_c:				1000.0
 * b_c:                         1000.0
 * lambda:                         1.0
 *
 * s1:				10
 * s2:				12.5
 *
 * alpha:  
 * beta:  
 * gamma: 
 * ...
 *
 * You can add more variables after the last one; 
 * IF YOU ADD MORE PARAMETERS, REMEMBER THE FOLLOWING:
 * 1) Increase the definition NUMBER_OF_PARAMETERS in read_parameter_file.h
 * 2) Add as many "fscanf(file_id, "%s %lf\n", header, &parameters[k] &parameters[k+1]);"
 *     as needed inside this function.
 * 3) modify "mountain_transformations.c" and "external_geometries.c" to read the new values correctly!
 *
 * simone.marras@gmail.com
 *******************************************************************/

int read_parameter_file(double *parameters, char *vertical_coords[])
{
  int error_flag; // error_flag==0 if the file does not exsist in the working directory
  
  unsigned int count;
  int i,j,k=0, file_nlines, line_cntr=0;
  char prob[12], prob_type[8], header[24];
  char *ptr_str;
  
  FILE *file_id;
  
  if( (file_id = fopen("mountain_parameters.inp", "r")) == NULL) {
    printf("Error opening file: 'mountain_parameters.inp'. %s\n", strerror(errno));
    error_flag = 0;
    exit (-1);
  }
  else{
    error_flag = 1;
    
    fscanf(file_id, "%s %s %s %s %s\n", header, header, header, header, header);
    fscanf(file_id, "%s %s %s %s %s %s\n", header, header, header, header, header, header);
    fscanf(file_id, "%s %s %s %s %s %s %s\n", header, header, header, header, header, header, header);
    
    //Read in the numerical entries from the file
    fscanf(file_id, "%s %lf %lf\n", header, &parameters[1], &parameters[2]);
    fscanf(file_id, "%s %lf %lf\n", header, &parameters[10], &parameters[11]);
    fscanf(file_id, "%s %lf\n",     header, &parameters[3]);
    fscanf(file_id, "%s %lf\n",     header, &parameters[4]); //hm or r
    fscanf(file_id, "%s %lf\n",     header, &parameters[5]); //a_c, or xc
    fscanf(file_id, "%s %lf\n",     header, &parameters[9]); //b_c, or yc
    fscanf(file_id, "%s %lf\n",     header, &parameters[6]); //lambda or yc
    fscanf(file_id, "%s %lf\n",     header, &parameters[7]); //s1 or theta0
    fscanf(file_id, "%s %lf\n",     header, &parameters[8]); //s2 or theta1
    fscanf(file_id, "%s %s\n",      header, vertical_coords[0]); //keyword for the vertical type of grid
    //END reading numerical values.
    
  }//End else.
  
  //Close file	
  fclose(file_id);
  
  printf(" #\n # MOUNTAIN: Entries from mountain_parameters_file\n");
  printf(" # Mountain Height:           %lf\n", parameters[4]);
  printf(" # Mountain 1/2 width in x:   %lf\n", parameters[5]);
  printf(" # Mountain 1/2 width in y:   %lf\n", parameters[9]);
  printf(" # s1:                        %lf\n", parameters[7]);
  printf(" # s2:                        %lf\n", parameters[8]);
  printf(" # Vertical coordinate::      %s\n", vertical_coords[0]);
  printf(" # END Mountain Entries from mountain_parameters_file\n");

  return error_flag;
}

/****************************************************************************************
 *
 * int CountNRows(char *input_file)
 *
 * This function counts the number of rows of a file and stops at the last numeric value.
 * If your file terminates with more than one blank line, this function does NOT 
 * count those additional lines.
 * simone.marras@gmail.com
 *
 *****************************************************************************************/
int CountNRows(char *input_file)
{
  int nrows;
  
  int  i, j, cntr, line_cntr;
  char line[256];
  
  FILE *file_ID;
  
  /********************************************************************************
   * Open file and read it throughout once:
   * FIRST opening of the file to read the number of elements for each category (GEOMETRY, GRID, etc.)
   * File *. 
   ********************************************************************************/
  if((file_ID = fopen(input_file, "r")) == NULL){
    printf(" ERROR in CountNRows ! The file %s could not be open\n", input_file);
    printf(" The program will EXIT now.\n\n");
    exit(1);
  }
  else{
   
    rewind(file_ID); //This line resets the pointer to the beginning of the file
    
    line_cntr = 0;
    while(fgets(line, sizeof(line), file_ID)!= NULL) {
      //zero, not Oh: This control helps to avoid to count any eventual empty line at the end.
      //if(line[0] < '0')
      //break;
      
      line_cntr++;
      
      // ELEMENTS
      nrows = line_cntr ;
      //End else
    }//End WHILE
  }
  fclose(file_ID);
  
  return nrows;
} //END CountNRows()


/*********************************************
 * load.c
 *
 * (function that reads a numeric file of matricial shape
 * e.g. 9 4 2
 *      2 3 5
 *      ... ... ...
 *
 * Similar to the Matlab function by the same name.)
 *
 * simone.marras@gmail.com
 **********************************************/

#define MAX_NODESPERDOMAIN 8

//Read a file that contains a matrix of floats (only up to 6 columns for now)
void load_input_float(char *file_name, float **input_domainNodes, int ncols)
{

  int count;
  int  i, j, k, cntr, line_cntr;
	
  FILE *file_ID;

  if((file_ID = fopen(file_name, "r")) == NULL)
    printf("The file %s could not be open\n", file_name);
  else{
    k=1;
    while(!feof(file_ID)){
      if(ncols == 1)
	fscanf(file_ID, "%f\n", &input_domainNodes[k][1]);
      else if(ncols == 2)
	fscanf(file_ID, "%f %f\n", &input_domainNodes[k][1], &input_domainNodes[k][2]);
      else if(ncols == 3)
	fscanf(file_ID, "%f %f %f\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3]);
      else if(ncols == 4)
	fscanf(file_ID, "%f %f %f %f\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4]);
      else if(ncols == 5)
	fscanf(file_ID, "%f %f %f %f %f\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4], &input_domainNodes[k][5]);
      else if(ncols == 6)
	fscanf(file_ID, "%f %f %f %f %f %f\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4], &input_domainNodes[k][5], &input_domainNodes[k][6]);
      else{
	printf(" !!! ERROR:\n !!! You are trying to read an input file with more columns than accepted.\n !!! SOLUTION: Open 'load_input.c' and add more ncols in the IF-ELSE section to add more\n\n");
	exit(1);
      }
			
      k=k+1;
    }
  }

  //Close file:
  fclose(file_ID);
		
		
  //free_matrix(input_domainNodes, 1,MAX_NODESPERDOMAIN, 1,2); // FREED in the calling function.

  return;
}


//Read a file that contains a matrix of integers (only up to 6 columns for now)
void load_input_int(char *file_name, int **input_domainNodes, int ncols)
{

  int count;
  int  i, j, k, cntr, line_cntr;
	
  FILE *file_ID;

  if((file_ID = fopen(file_name, "r")) == NULL)
    printf("The file %s could not be open\n", file_name);
  else{
    k=1;
    while(!feof(file_ID)){
      if(ncols == 1)
	fscanf(file_ID, "%d\n", &input_domainNodes[k][1]);
      else if(ncols == 2)
	fscanf(file_ID, "%d %d\n", &input_domainNodes[k][1], &input_domainNodes[k][2]);
      else if(ncols == 3)
	fscanf(file_ID, "%d %d %d\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3]);
      else if(ncols == 4)
	fscanf(file_ID, "%d %d %d %d\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4]);
      else if(ncols == 5)
	fscanf(file_ID, "%d %d %d %d %d\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4], &input_domainNodes[k][5]);
      else if(ncols == 6)
	fscanf(file_ID, "%d %d %d %d %d %d\n", &input_domainNodes[k][1], &input_domainNodes[k][2], &input_domainNodes[k][3], &input_domainNodes[k][4], &input_domainNodes[k][5], &input_domainNodes[k][6]);
      else{
	printf(" !!! ERROR:\n !!! You are trying to read an input file with more columns than accepted.\n !!! SOLUTION: Open 'load_input.c' and add more ncols in the IF-ELSE section to add more\n\n");
	exit(1);
      }
			
      k=k+1;
    }
  }

  //Close file:
  fclose(file_ID);
		
		
  //free_imatrix(input_domainNodes, 1,MAX_NODESPERDOMAIN, 1,2); // FREED in the calling function.

  return;
}


//Same thing but for 1D arrays:
void load_input_fvect(char *file_name, float *input_vect)
{
  
  int count;
  int  i, j, k, cntr, line_cntr;
  
  FILE *file_ID;
  
  if((file_ID = fopen(file_name, "r")) == NULL)
    printf("The file %s could not be open\n", file_name);
  else{
    k=1;
    while(!feof(file_ID)){
      fscanf(file_ID, "%f\n", &input_vect[k]);
    }
    
    k=k+1;
  }
  
  
  //Close file:
  fclose(file_ID);
  
  return;
}



void load_input_ivect(char *file_name, int *input_vect)
{

  int count;
  int  i, j, k, cntr, line_cntr;
  
  FILE *file_ID;
  
  if((file_ID = fopen(file_name, "r")) == NULL)
    printf("The file %s could not be open\n", file_name);
  else{
    k=1;
    while(!feof(file_ID)){
      fscanf(file_ID, "%d\n", &input_vect[k]);
    }
    
    k=k+1;
  }
  
  
  //Close file:
  fclose(file_ID);
  
  return;
}


