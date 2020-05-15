/***************************************************************************
 *
 * topo_user_function.c
 *
 * This function builds a surface with an analytic function.
 * The function to be used is defined by the number (itask)
 * given in input by the user if the TOPOGRAPHY entry is "FUNCTION"
 *
 * The number is read in READ_INPUT.c and stored in problem[5] 
 *
 ***************************************************************************/

#include<string.h>
#include<stdio.h>
#include<math.h>

#include "global_vars.h"

double topo_user_3Dfunction(char *fun_type, double x, double y, float a, float b, float height, float coeff)
{
  double f;
  double a_c,hm;
  double x2, x3, x4;
  
  if( !strcmp(fun_type, "0") )
    f = 0.0;

  else if( !strcmp(fun_type, "1") )
    f = height/pow(1.0 + (x-xc)*(x-xc)/a/a +  (y-yc)*(y-yc)/b/b, coeff);

  else if( !strcmp(fun_type, "2") )
    f = height*sin(x-xc)*cos(y-yc)/(sqrt(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))) + 1.0)*10.0;
  
  else if( !strcmp(fun_type, "3") )
    {
      //Compute the mountain coordinates:
      a_c = xmax*0.65; //so that we don't worry on the input values.
      hm  = xmax*0.35;
      f   = -(a_c + hm*sin(PI*(0.5 + 2.0/xmax * x)) - zmax); //Ok
    }
  else if(!strcmp(fun_type, "4") )
    {
      x2 = x*x;
      x3 = x*x2;
      x4 = x*x3;
      f=  0.594689181*(0.298222773*sqrt(x) - 0.127125232*x - 0.357907906*x2 + 0.291984971*x3 - 0.105174606*x4);
      
    }
  else
    printf(" ERROR in function definition number: only FUNCTION 1 and 2 are programmed\n");
    
  return f;
}
