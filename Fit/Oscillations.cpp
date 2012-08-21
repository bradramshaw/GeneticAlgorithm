#include "StdAfx.h"
#include "Oscillations.h"

#define PI 3.1415926535897932384626433832795028841971
namespace Oscillations
{
 double Oscillations::oscillations(double x, double F, double dF, double dF2, double H)
	{
		double y = sin(2*PI*(F+dF*cos(x)+dF2*cos(2*x))/H);
		return y;
	}


}