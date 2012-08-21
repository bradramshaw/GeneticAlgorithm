#pragma once
#include "Envelope.h"
#include "OscillationFunction.h"
#include "Integrator.h"
#include "Parameters.h"
#include "FitThreadFast.h"

//#include <boost\random\mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>

class RandomFitThreadedFast
{
	public:
	
		FitThreadFast ** p_threads;
		Parameters::fitParameters min_parameters;
			
		double ** data_set;
	
		int length;		
		int iterations;
		int n_threads;
	
		RandomFitThreadedFast(void);
		RandomFitThreadedFast(double ** _data_set, int _length, int _n_threads, int _iterations);
		~RandomFitThreadedFast(void);


		double getMinimum();

		void setIterations(int _iterations);
		void setThreadNumber(int _threads);

		void outputParams();
};

