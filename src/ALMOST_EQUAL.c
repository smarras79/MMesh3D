#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "ALMOST_EQUAL.h"

#define TOL 2*DBL_EPSILON

bool AlmostEqual(double a, double b)
{
    bool isalmostEqaul;

    if (a == 0.0 || b == 0.0)
	{
	    if (fabs(a - b) <= TOL)
		isalmostEqaul = true;
	    else
		isalmostEqaul = false;	    
	}
    else
	{
	    if (fabs(a - b) <= TOL*fabs(a) && fabs(a - b) <= TOL*fabs(b))
		isalmostEqaul = true;
	    else
		isalmostEqaul = false;
	}
        
    return isalmostEqaul;
}
