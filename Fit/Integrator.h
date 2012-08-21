#pragma once
class Integrator
{
private:


public:
	

	static const double weights[];
	static const double nodes[];
	
	Integrator();
	~Integrator(void);
	//double integrate(void* function, double (*functionMethod)(void*,double), double lower_limit, double upper_limit, int sample_points);
//	static double integrateSimpsons(void* function, double (*functionMethod)(void*,double), double lower_limit, double upper_limit, int sample_points);

	static double integrateLegendre( void* function, double (*functionMethod)( void*,int, double, double), double H, double theta, double lower_limit, double upper_limit);
	
	
	
	
	/*static double integrateLegendre2( double H, double theta);

	static const double coeffs[];
	static const double n_coeffs;
	double cos_vals[n_nodes];
	double cos2_vals[n_nodes];
	double F,dF,dF2,phi;*/
};

