/*******************************************************************
 * minmaxval.c
 * 
 * Function that finds the maximum and minimum values inside an array
 *
 * simone.marras@gmail.com
 *******************************************************************/

#include<stdio.h>

float maxval1Darray( float *ARRAY, int nrows, int index)
{
  int i;
  float maxvalue, minvalue;
  
  maxvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++) {
    
    if(ARRAY[i] > maxvalue ){
      maxvalue = ARRAY[i];
    }
  }
  
  return maxvalue;
}

float minval1Darray( float *ARRAY, int nrows, int index)
{
  
  int i;
  float maxvalue, minvalue;
  
  minvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++)
    {
      if (ARRAY[i] < minvalue )
	{
	  minvalue = ARRAY[i];
	}    
    }
  
  return minvalue;
}

float maxval2Dfloat( float **ARRAY, int nrows, int column, int index)
{
  int i;
  float maxvalue, minvalue;
  
  //Find the horizontal domain limits (in -x and +x):
  maxvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++) {
    
    //Computes the maximum of z-coordinates directly in this loop:
    if(ARRAY[i][column] > maxvalue ){
      maxvalue = ARRAY[i][column];
    }
  }
  
  return maxvalue;
}

float minval2Dfloat( float **ARRAY, int nrows, int column, int index)
{
  
  int i;
  float maxvalue, minvalue;
  
  //Find the horizontal domain limits (in -x and +x):
  minvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++)
    {
      if (ARRAY[i][column] < minvalue )
	{
	  minvalue = ARRAY[i][column];
	}    
    }
  
  return minvalue;
}

double maxval2Ddouble( double **ARRAY, int nrows, int column, int index)
{
  int i;
  double maxvalue, minvalue;
  printf(" nrows=%d\n", nrows);
  //Find the horizontal domain limits (in -x and +x):
  maxvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++) {
  
    //Computes the maximum of z-coordinates directly in this loop:
    if(ARRAY[i][column] > maxvalue ){
      maxvalue = ARRAY[i][column];
    }
  }
 
  return maxvalue;
}

double minval2Ddouble( double **ARRAY, int nrows, int column, int index)
{
  
  int i;
  double maxvalue, minvalue;
  
  //Find the horizontal domain limits (in -x and +x):
  minvalue = 0.0;
  for(i = index+1; i < nrows-1+index; i++)
    {
      if (ARRAY[i][column] < minvalue )
	{
	  minvalue = ARRAY[i][column];
	}    
    }
  
  return minvalue;
}

