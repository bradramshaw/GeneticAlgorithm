#pragma once
#define n_nodes 60
#define n_nodesp5 30
#define OSC_CONST 14.693185881921599
class IntegratorFast
{
public:

	static const double weights[];
	static const double nodes[];
	static const double coeffs[];
	static const int n_coeffs;
	double cos_vals[n_nodes];
	double cos2_vals[n_nodes];
	double F,dF,dF2,phi, A,T,Td,m;

	IntegratorFast(double _F,double _dF,double _dF2, double _phi,double _A,double _T,double _Td, double _m);
	~IntegratorFast(void);


	double integrateLegendre( double H, double theta, double lower_limit, double upper_limit);
	double getYValue(int x, double H, double theta);
	double besselJ_C(double _x);
	void changeParams(double _F,double _dF,double _dF2, double _phi,double _A,double _T,double _Td, double _m);
};

