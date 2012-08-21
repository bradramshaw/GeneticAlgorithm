#include "StdAfx.h"
#include "Function.h"

Function::Function(void)
{

}

Function::~Function(void)
{

}
double Function::getYValue(double x)
{
	double y = x;
	return y;
}


/*Function::Function(double * x_vals, int points)
{
	Function::y_calculated = 0;
	Function::x_values = x_vals;
	Function::y_values = new double[points];
	Function::number_of_points = points;
}


Function::~Function(void)
{
	delete[] Function::x_values;
	delete[] Function::y_values;
}

void Function::calculateYValues(int number_of_points)
{
	for(int i = 0; i < number_of_points; i++)
	{
		Function::y_values[i] = Function::x_values[i];
	}

}
double * Function::getYValues()
{
	if( y_calculated == 0)
	{
		Function::calculateYValues(number_of_points);
		y_calculated = 1;
	}
	return Function::y_values;
}

double * Function::getXValues()
{
	return Function::x_values;
}

int Function::getNumberOfPoints()
{
	return Function::number_of_points;
}*/