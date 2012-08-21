#include "StdAfx.h"
#include "OscillationFunction.h"


#define PI 3.1415926535897932384626433832795028841971
#define PI2 6.28318530717958647692528676656
#define besselj boost::math::cyl_bessel_j<int, double>
#define KF_CONV 0.0644735999999999
#define tan_0 0
#define cos_0 1

const double OscillationFunction::nodes[] ={-0.999210123227436,-0.995840525118838,-0.989787895222222,-0.981067201752598,-0.969701788765053,-0.955722255839996,-0.939166276116423,-0.920078476177628,-0.898510310810046,-0.874519922646899,-0.84817198478593,-0.819537526162146,-0.788693739932264,-0.7557237753065849,-0.7207165133557299,-0.683766327381355,-0.644972828489477,-0.60444059704851,-0.562278900753944,-0.518601400058569,-0.473525841761707,-0.427173741583078,-0.379670056576798,-0.331142848268447,-0.281722937423261,-0.231543551376029,-0.180739964873425,-0.129449135396944,-0.07780933394953649,-0.0259597723012478,0.0259597723012475,0.07780933394953619,0.129449135396945,0.180739964873425,0.231543551376029,0.281722937423261,0.331142848268448,0.379670056576797,0.427173741583078,0.473525841761707,0.518601400058569,0.562278900753944,0.60444059704851,0.644972828489477,0.683766327381355,0.7207165133557299,0.7557237753065849,0.788693739932264,0.819537526162146,0.84817198478593,0.874519922646899,0.898510310810046,0.920078476177627,0.939166276116424,0.955722255839996,0.969701788765053,0.981067201752598,0.989787895222222,0.995840525118838,0.999210123227436};
const double OscillationFunction::coeffs[] = {1.00000000001502,-0.0000000003971037387806575,-0.2500000025918557,0.00000004370805323228311,0.01562481926022373,0.0000003746964069017977,-0.0004344981963593569,0.0000003923846648059681,0.000006552041596216036,0.000000097577427433624,-0.0000000985582498286686,0.000000007259332734370126,-0.000000000815801230627031,0.0000000001695776376652846,-0.00000000001859116735850654,0.000000000001055434172075196,-0.00000000000003087777334068834,0.0000000000000003739549721015043};
const int OscillationFunction::n_coeffs = 18;

OscillationFunction::OscillationFunction(double _F,double _dF,double _dF2, double _phi) : Function()
{
	F = _F;
	dF = _dF;
	dF2 = _dF2;
	phi = _phi;
	for(int i =0; i<n_nodes; i++)
	{
		OscillationFunction::cos_vals[i] = cos(PI*nodes[i]);
		OscillationFunction::cos2_vals[i] = cos(2*PI*nodes[i]);
	}
	
}

double OscillationFunction::getYValue(int x, double H, double theta)
{

	
	double kf_sqF_tan0 = KF_CONV * sqrt(F) * tan_0;
	return sin(PI2*((F+dF*cos_vals[x]*besselJ_C(kf_sqF_tan0) + dF2*cos2_vals[x]*besselJ_C(2*kf_sqF_tan0))/(H*cos_0)-phi));
	//return sin(2*PI*((F+dF*cos_vals[x]*besselj(0,KF_CONV * sqrt(F) * tan(theta)) + dF2*cos2_vals[x]*besselj(0,2*KF_CONV * sqrt(F) * tan(theta)))/(H*cos(theta))-phi));
	//double y = sin(2*PI*((F+dF*cos_vals[x]*besselJ_C(KF_CONV * sqrt(F) * tan(theta)) + dF2*cos2_vals[x]*besselJ_C(2*KF_CONV * sqrt(F) * tan(theta)))/(H*cos(theta))-phi));
	
}

inline double OscillationFunction::besselJ_C(double _x)
{
		
	//double y =1.00000000001502-0.000000000397103738780658*_x-0.250000002591856*_x*_x+0.0000000437080532322831*_x*_x*_x+0.0156248192602237*_x*_x*_x*_x+0.000000374696406901798*_x*_x*_x*_x*_x-0.000434498196359357*_x*_x*_x*_x*_x*_x+0.000000392384664805968*_x*_x*_x*_x*_x*_x*_x+0.00000655204159621604*_x*_x*_x*_x*_x*_x*_x*_x+0.000000097577427433624*_x*_x*_x*_x*_x*_x*_x*_x*_x-0.0000000985582498286686*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x+0.00000000725933273437013*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x-0.000000000815801230627031*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x+0.000000000169577637665285*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x-0.0000000000185911673585065*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x+0.0000000000010554341720752*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x-0.0000000000000308777733406883*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x+0.000000000000000373954972101504*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x*_x;
	//return y;

	double total = 0;
	double x_pow = 1;
	int i = 0;
	
	do
	{
		total += coeffs[i]*x_pow;
		x_pow *= _x;
		i++;
	}
	while(i < n_coeffs);

	return total;

	/*double total = 0;
	double x_pow = 1;
	const double * coeff_add = &coeffs[0];
	
	do
	{
		total += (*coeff_add)*x_pow;
		x_pow *= _x;
		coeff_add++;
	}
	while(coeff_add < &coeffs[n_coeffs-1]);

	return total;*/
}

double OscillationFunction::wrapper_getYValue( void* pointer, int x, double H, double theta)
{
	
	 OscillationFunction * myFunction = (OscillationFunction*) pointer;
	return myFunction->getYValue(x, H, theta);
	
}

void OscillationFunction::changeParams(double _F,double _dF,double _dF2, double _phi)
{
	F = _F;
	dF = _dF;
	dF2 = _dF2;
	phi = _phi;
}

