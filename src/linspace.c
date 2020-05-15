/* linspace.c
 *
 * Vfloat  = flinspace(vmin, vmax, n);
 * Vint    = ilinspace(vmin, vmax, n);
 * Vdouble = dlinspace(vmin, vmax, n);
 * as linspace in Matlab and Octave: it defines an array of values from vmin to vmax with n values.
 *
 * WARNING: THESE FUNCTIONS ARE WRITTEN TO ACCEPT 1-indexed arrays (instead of the 0-indexed as in standard C arrays)!!!
 */

#include<stdio.h>
#include"nrutil.h"

//Floats
void flinspace(float vmin, float vmax, int n, float *VECT, char *force_flag)
{
	//!!! VECT must be allocated starting from 1 for this function to work properely
	
	int i;
	float dv;

	VECT[1] = vmin;
	for (i=2; i<=n; i++){
		dv = (vmax - vmin)/(n-1);
		VECT[i] = VECT[i-1] + dv;
	}
	
	 // Force the last point from linspace to coincide with the extreme that was chosen:
	if( !strcmp(force_flag,"f") || !strcmp(force_flag,"force") ){
		if(VECT[n] != vmax){
			VECT[n] = vmax;
		}
	}
	//End forcing coincidence of values.)

	//printf(" If not done already, remember to deallocate this VECT within the calling function\n");

return;
}

//Integers
void ilinspace(int vmin, int vmax, int n, int *VECT)
{
	//!!! VECT must be allocated starting from 1 for this function to work properely
	
	int i;
	float dv;
	
	VECT[1] = vmin;
	for (i=2; i<=n; i++){
		dv = (vmax - vmin)/(n-1);
		dv = (int)dv;
		VECT[i] = VECT[i-1] + dv;
	}

	//printf(" If not done already, remember to deallocate this VECT within the calling function\n");

return;
}

//Double
void dlinspace(double vmin, double vmax, int n, double *VECT, char *force_flag)
{
	//!!! VECT must be allocated starting from 1 for this function to work properely
	
	int i;
	double dv;

	VECT[1] = vmin;
	for (i=2; i<=n; i++){
		dv = (vmax - vmin)/(n-1);
		VECT[i] = VECT[i-1] + dv;
	}
	
	 // Force the last point from linspace to coincide with the extreme that was chosen:
	if( !strcmp(force_flag,"f") || !strcmp(force_flag,"force") ){
		if(VECT[n] != vmax){
			VECT[n] = vmax;
		}
	}
	//End forcing coincidence of values.)

	//printf(" If not done already, remember to deallocate this VECT within the calling function\n");

return;
}


void linspace_test(float vmin, float vmax, int n, float *VECT)
{
	//!!! VECT must be allocated starting from 1 for this function to work properely
	
	int i;
	double dv;
printf(" vmin: %f   vmax:%f\n",vmin,vmax);
	VECT[1] = vmin;
	for (i=2; i<=n; i++){
		dv = (vmax - vmin)/(n-1);
		VECT[i] = VECT[i-1] + dv;
		printf("dv=%f, VECT[%d]: %f\n",dv, i,VECT[i]);
	}

	//printf(" If not done already, remember to deallocate this VECT within the calling function\n");

return;
}
