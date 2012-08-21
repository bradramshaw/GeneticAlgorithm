#pragma once
#include "Function.h"

#define n_nodes 60

class OscillationFunction: public Function
{
public:
	static const double nodes[];
	static const double coeffs[];
	static const int n_coeffs;
	double cos_vals[n_nodes];
	double cos2_vals[n_nodes];
	double F,dF,dF2,phi;
	OscillationFunction(double _F,double _dF,double _dF2, double _phi);
	OscillationFunction();
	void changeParams(double _F,double _dF,double _dF2, double _phi);
	double getYValue(int x, double H, double theta);
	static double wrapper_getYValue( void* pointer, int x, double H, double theta);
	static double besselJ_C(double _x);
};
