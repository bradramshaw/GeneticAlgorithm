#include "StdAfx.h"
#include "Envelope.h"
#define OSC_CONST 14.693185881921599

Envelope::Envelope(double _A,double _T,double _Td, double _m) : Function()
{
	A = _A;
	T = _T;
	Td = _Td;
	m = _m;
	
}

double Envelope::getYValue(double H, double theta)
{
	double y = A*sqrt(H)*exp(-OSC_CONST*Td*m/(H*cos(theta)))*(OSC_CONST*T*m/(H*cos(theta)))/sinh(OSC_CONST*T*m/(H*cos(theta)));
	return y;
}

void Envelope::changeParams(double _A,double _T,double _Td, double _m)
{
	A = _A;
	T = _T;
	Td = _Td;
	m = _m;
	
}