#pragma once

class Function
{
public:
	Function(void);
	~Function(void);
	static double getYValue(double x);
};

/*class Function
{
public:
	Function(double * x_vals, int number_of_points);
	~Function(void);
	virtual double * getYValues();
	double * getXValues();
	int Function::getNumberOfPoints();

protected: 
	int y_calculated;
	double * x_values;
	double * y_values;
	int number_of_points;
	virtual void calculateYValues(int number_of_points);
};*/

