#pragma once
#include "Envelope.h"
#include "OscillationFunction.h"
#include "Integrator.h"
#include "Parameters.h"
#include "FitThread.h"

//#include <boost\random\mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>

class RandomFitThreaded
{
	public:
	
		FitThread ** p_threads;
		Parameters::fitParameters min_parameters;
			
		double ** data_set;
	
		int length;		
		int iterations;
		int n_threads;
	
		RandomFitThreaded(void);
		RandomFitThreaded(double ** _data_set, int _length, int _n_threads, int _iterations);
		~RandomFitThreaded(void);


		double getMinimum();

		void setIterations(int _iterations);
		void setThreadNumber(int _threads);

		void outputParams();
};

