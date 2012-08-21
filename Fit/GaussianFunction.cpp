#include "StdAfx.h"
#include "GaussianFunction.h"

double GaussianFunction::getYValue(double x)
{
	double y = exp(-pow(x,2));
	return y;
}


/*
 GaussianFunction::GaussianFunction(double * x_vals, int points): Function(x_vals, points)
{
	
}
 
void GaussianFunction::calculateYValues(int number_of_points)
{
	for(int i = 0; i < number_of_points; i++)
		{
			y_values[i] = exp(-pow(
			
		}
}

double * GaussianFunction::getYValues()
{
	if( y_calculated == 0)
	{
		GaussianFunction::calculateYValues(number_of_points);
		y_calculated = 1;
	}
	return Function::y_values;
}*/
