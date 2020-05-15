/* File: read_geogrid.c

   Sample subroutine to read an array from the geogrid binary format.

   Notes: Depending on the compiler and compiler flags, the name of 
   the read_geogrid() routine may need to be adjusted with respect
   to the number of trailing underscores when calling from Fortran.

   Michael G. Duda, NCAR/MMM
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define SEPARATOR "= "  //Used by READ_INDEX()

//#ifdef _UNDERSCORE
//#define read_geogrid read_geogrid_
//#endif
//#ifdef __DOUBLEUNDERSCORE
//#define read_geogrid read_geogrid__
//#endif

//#define BIG_ENDIAN    0
//#define LITTLE_ENDIAN 1

#define big_endian    0
#define little_endian 1

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
      int * status)
{
   int i, ival, cnt, narray;
   int A2, B2;
   int A3, B3, C3;
   int A4, B4, C4, D4;
   unsigned char * c;
   char local_fname[1024];
   FILE * bfile;

   *status = 0;

   narray = (*nx) * (*ny) * (*nz);

   /* Make a null-terminated local copy of the filename */
   strncpy(local_fname,fname,*len);
   local_fname[*len]='\0';
   
   /* Attempt to open file for reading */
   if (!(bfile = fopen(local_fname,"rb")))
   {
     printf(" Geo file %s could not be open\n", local_fname);
      *status = 1;
      return 1;
   }

   /* Allocate memory to hold bytes from file and read data */
   c = (unsigned char *)malloc(sizeof(unsigned char)*(*wordsize) * narray);
   cnt = fread((void *)c, sizeof(unsigned char), narray*(*wordsize), bfile);
 
   fclose(bfile);

   if (cnt == 0) 
   {
      *status = 1;
      return 1;
   }

   /*
      Set up byte offsets for each wordsize depending on byte order.
      A, B, C, D give the offsets of the LSB through MSB (i.e., for 
      word ABCD, A=MSB, D=LSB) in the array from the beginning of a word 
   */
   if (*endian == big_endian) {
      A2 = 0; B2 = 1;
      A3 = 0; B3 = 1; C3 = 2;
      A4 = 0; B4 = 1; C4 = 2; D4 = 3;
   }
   else {
      B2 = 0; A2 = 1;
      C3 = 0; B3 = 1; A3 = 2;
      D4 = 0; C4 = 1; B4 = 2; A4 = 3;
   }

   /* Convert words from native byte order */
   switch(*wordsize) {
      case 1:
         for(i=0; i<narray; i++)
         {
            ival = (int)(c[i]);      
            if ((*isigned) && (ival > (1 << 7))) ival -= (1 << 8);
            rarray[i] = (float)ival;
         }
         break;

      case 2:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[2*i+A2]<<8) | (c[2*i+B2]));      
            if ((*isigned) && (ival > (1 << 15))) ival -= (1 << 16);
            rarray[i] = (float)ival;
	    //printf("%d %f \n",i,rarray[i]); 
         }
         break;

      case 3:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[3*i+A3]<<16) | (c[3*i+B3]<<8) | c[3*i+C3]);      
            if ((*isigned) * (ival > (1 << 23))) ival -= (1 << 24);
            rarray[i] = (float)ival;
         }
         break;

      case 4:
         for(i=0; i<narray; i++)
         {
            ival = (int)((c[4*i+A4]<<24) | (c[4*i+B4]<<16) | (c[4*i+C4]<<8) | c[4*i+D4]);      
            if ((*isigned) && (ival > (1 << 31))) ival -= (1 << 32);
            rarray[i] = (float)ival;
	    //printf("%d %f \n",i,rarray[i]); 
         }
         break;
   }

   free(c);

   /* Scale real-valued array by scalefactor */
   if (*scalefactor != 1.0)
     {
       for (i=0; i<narray; i++)
	 rarray[i] = rarray[i] * (*scalefactor);
     }
   
   //Print rarray:
   //for (i=0; i<narray; i++)   
   // printf("%d %f \n",i,rarray[i]); 
  
   return 0;
}

/*****************************************************************************
 * Function to read the "index" description file of the topography files:
 *
 * simone.marras@bgmail.com
 * May 22, 2008
 *
 *****************************************************************************/

void READ_INDEX(char *input_inp, char **topo_type, char **topo_sign, char **topo_projection, char **topo_units, char **topo_description, double *dx, double *dy, double *known_x, double *known_y, double *known_lat, double *known_lon, int *wordsize, int *tile_x, int *tile_y, int *tile_z, int *tile_bdr)
{
  
  int line_cntr, cntr;
  int type_count, signed_count, proj_count, dx_count, dy_count;
  int knownx_count, knowny_count, knownz_count, knownlat_count, knownlon_count;
  int wordsz_count, tilex_count, tiley_count, tilez_count, tilebdr_count, units_count, descr_count;
  
  char line[128];
  char *ptr_str;
  
  FILE *input_file_id;
  
  /***************************************************************************************
   *1st OPENING Of the FILE:
   ***************************************************************************************/
  printf("\nThe user input_file is '%s'\n...searching for the file...\n\n", input_inp);
  
  if((input_file_id = fopen(input_inp, "r")) == NULL){
    printf("The input file could not be open.\nVerify that it exists in the current directory\nThe program will exit now\n\n");
    return;
  }
  else{
    printf("!File '%s' opened successfully\n\n", input_inp);
    
  }
  
  rewind(input_file_id); //This line resets the pointer to the beginning of the file
  line_cntr = 0;
  while(fgets(line, sizeof(line), input_file_id)!= NULL) {
    
    line_cntr++;
    
    // Part to sepatrate word by word each line of the file:
    ptr_str = strtok(line, SEPARATOR);	//This line breaks the original line into
    //as many "words" there are in that line.
    //SEPARATOR tells "strtok" how the strings are separted.
    
    /********************************************************************************************/	
    //In the new few lines we make a control over the different heading appearing in the index file:
    // type, signed, dx, etc.
    //This is necessary to count the numeber of inputs that correspond to each category.
    //This solution makes this code very versatile if the number of entry in each category is changed by the user!
    
    //0) type=
    if(strncmp(ptr_str,"type", 4) == 0){		
      type_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], type_count);
    }
    //1) signed=:
    if(strncmp(ptr_str,"signed", 6) == 0){		
      signed_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], signed_count);
    }
    //2) projection=:
    if(strncmp(ptr_str,"projection", 10) == 0){		
      proj_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], proj_count);
    }
    //3) dx=:
    if(strncmp(ptr_str,"dx", 2) == 0){		
      dx_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], dx_count);
    }
    //4) dy=:
    if(strncmp(ptr_str,"dy", 2) == 0){		
      dy_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], dy_count);
    }	
    //5) known_x:
    if(strncmp(ptr_str,"known_x", 7) == 0){		
      knownx_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], knownx_count);
    }
    //6) known_y:
    if(strncmp(ptr_str,"known_y", 7) == 0){		
      knowny_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], knowny_count);
    }
    //7) know_lat:
    if(strncmp(ptr_str,"known_lat", 9) == 0){		
      knownlat_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], knownlat_count);
    }
    //8) know_lon:
    if(strncmp(ptr_str,"known_lon", 9) == 0){		
      knownlon_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], knownlon_count);
    }
    //9) wordsize:
    if(strncmp(ptr_str,"wordsize", 8) == 0){		
      wordsz_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], wordsz_count);
    }
    //10) tile_x:
    if(strncmp(ptr_str,"tile_x", 6) == 0){		
      tilex_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], tilex_count);
    }
    //11) tile_y:
    if(strncmp(ptr_str,"tile_y", 6) == 0){		
      tiley_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], tiley_count);
    }
    //12) tile_z:
    if(strncmp(ptr_str,"tile_z", 6) == 0){		
      tilez_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], tilez_count);
    }
    //13) tile_bdr:
    if(strncmp(ptr_str,"tile_bdr", 6) == 0){		
      tilebdr_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], tilebdr_count);
    }
    //14) units:
    if(strncmp(ptr_str,"units", 5) == 0){		
      units_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], units_count);
    }
    //15) description:
    if(strncmp(ptr_str,"description", 8) == 0){		
      descr_count = line_cntr;
      //printf("Category %s starts at line %d\n", &ptr_str[0], descr_count);
    }
    
  }//End WHILE
  
  //We NEED to close the file so that the counters of the category elements are not rewritten again.
  fclose(input_file_id);
  
  /***************************************************************************************************/
  /*2ns OPENING Of the FILE:
    /***************************************************************************************************/
  printf("\nThe user input_file is '%s'\n...searching for the file...\n\n", input_inp);
  
  if((input_file_id = fopen(input_inp, "r")) == NULL){
    printf("The input file could not be open.\nVerify that it exists in the current directory\nThe program will exit now\n\n");
    return;
  }
  else{
    printf("!File '%s' re-opened successfully for reading\n\n", input_inp);
    
  }
  
  rewind(input_file_id); //This line resets the pointer to the beginning of the file
  line_cntr = 0;
  while(fgets(line, sizeof(line), input_file_id)!= NULL) {
    
    line_cntr++;
    
    // Part to sepatrate word by word each line of the file:
    ptr_str = strtok(line, " ");	//This line breaks the original line into
    //as many "words" there are in that line.
    //SEPARATOR tells "strtok" how the strings are separted.
    
    // Part to read the actual values of the variables:
    
    //0) type:
    if(line_cntr == type_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  strcpy(*topo_type, ptr_str);
	  //printf("Type = %s\n", *topo_type);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (TYPE)
    //1) signed:
    if(line_cntr == signed_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  strcpy(*topo_sign, ptr_str);
	  //printf("Signed = %s\n", *topo_sign);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (SIGNED)
    //2) projection:
    if(line_cntr == proj_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  strcpy(*topo_projection, ptr_str);
	  //printf("Projection = %s\n", *topo_projection);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (PROJECTION)
    //3) dx:
    if(line_cntr == dx_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *dx = atof(ptr_str);
	  //printf("dx = %f\n", *dx);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (DX)
    //4) dy:
    if(line_cntr == dy_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *dy = atof(ptr_str);
	  //printf("dy = %f\n", *dy);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (DY)
    //5) knownx:
    if(line_cntr == knownx_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *known_x = atof(ptr_str);
	  //printf("known_x = %.2f\n", *known_x);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (KNOWN_X)
    //6) knowny:
    if(line_cntr == knowny_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *known_y = atof(ptr_str);
	  //printf("known_y = %.2f\n\n", *known_y);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (KNOWN_Y)
    //7) known_lat:
    if(line_cntr == knownlat_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *known_lat = atof(ptr_str);
	  //printf("known_lat = %f\n", *known_lat);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (KNOWN_LAT)
    //8) known_lon:
    if(line_cntr == knownlon_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *known_lon = atof(ptr_str);
	  //printf("known_lon = %f\n", *known_lon);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (KNOWN_LON)
    //9) wordsize:
    if(line_cntr == wordsz_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *wordsize = atof(ptr_str);
	  //printf("wordsize = %d\n", *wordsize);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (WORDSIZE)
    //10) tile_x:
    if(line_cntr == tilex_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *tile_x = atof(ptr_str);
	  //printf("tile_x = %d\n", *tile_x);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (TILE_X)
    //11) tile_Y:
    if(line_cntr == tiley_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *tile_y = atof(ptr_str);
	  //printf("tile_y = %d\n", *tile_y);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (TILE_Y)
    //12) tile_z:
    if(line_cntr == tilez_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *tile_z = atof(ptr_str);
	  //printf("tile_z = %d\n", *tile_z);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (TILE_Z)
    //13) tile_bdr:
    if(line_cntr == tilebdr_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  *tile_bdr = atof(ptr_str);
	  //printf("tile_bdr = %d\n", *tile_bdr);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (TILE_BDR)
    //14) units:
    if(line_cntr == units_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  strcpy(*topo_units,  ptr_str);
	  //printf("Units = %s\n", *topo_units);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (UNITS)
    //15) Description:
    if(line_cntr == descr_count){
      cntr = 0;
      while(ptr_str != NULL){
	switch(cntr){
	case 0 :
	  //printf("Variable name: %s\n", ptr_str);					
	  //printf("counter = %d\n", cntr);
	  break;
	case 1 :
	  strcpy(*topo_description, ptr_str);
	  //printf("Description = %s\n\n", *topo_description);
	  break;
	}
	
	ptr_str = strtok(NULL, SEPARATOR);
	cntr = cntr + 1;
      }//End WHILE
    }//End IF (DESCRIPTION)
    
    /***************************************************************************************************
     *CLOSING OF THE FILE after reading and assigning values to the variables:
     ***************************************************************************************************/
  }
  fclose(input_file_id);
  
  //Print the values all here: (comment them if you dont want them printed on screen:
  printf("####################################################\n");
  printf("# CONTENT of '%s' ::", input_inp);
  printf("\n# Topo_type = %s", *topo_type);
  printf("# Signed = %s", *topo_sign);
  printf("# Projection = %s", *topo_projection);
  printf("# dx = %f\n", *dx);
  printf("# dy = %f\n", *dy);
  printf("# known_x = %f\n", *known_x);
  printf("# known_y = %f\n", *known_y);
  printf("# known_lat = %f\n", *known_lat);
  printf("# known_lon = %f\n", *known_lon);
  printf("# wordsize = %d\n", *wordsize);
  printf("# Tile_x = %d\n", *tile_x);
  printf("# Tile_y = %d\n", *tile_y);
  printf("# Tile_z = %d\n", *tile_z);
  printf("# Tile_bdr = %d\n", *tile_bdr);
  printf("# Units = %s\n", *topo_units);
  printf("# Description = %s\n", *topo_description);
  printf("####################################################\n");
  printf("\n");
  
  printf("Reading of file '%s' has now completed\n", input_inp);
  
  
  return;
  
}//END READ_INDEX.c
