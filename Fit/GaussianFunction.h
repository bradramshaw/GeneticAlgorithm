#pragma once
#include "Function.h"



class GaussianFunction: public Function
{
public:
	static double getYValue(double x);
};



/*class GaussianFunction: public Function
{
public:
	GaussianFunction(double * x_vals, int points);
	virtual void calculateYValues(int number_of_points);
	virtual double * getYValues();
};*/