#pragma once
#include "function.h"
class Envelope : public Function
{
public:
	Envelope(double _A,double _T,double _Td,double _m);
	Envelope();
	double A,T,Td,m;
	void changeParams(double _A,double _T,double _Td,double _m);
	double getYValue(double H, double theta);
};

