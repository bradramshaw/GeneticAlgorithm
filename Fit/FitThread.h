#pragma once
#include "Envelope.h"
#include "OscillationFunction.h"
#include "Integrator.h"
#include "Parameters.h"

class FitThread
{
public:
	struct threadContents
	{
		HANDLE hEvent;
		int n_threads;
		int thread_id;
		int iterations;
		OscillationFunction *myOscFun;
		Envelope *myEnvel;
		Integrator *myInt;
		Parameters::fitParameters myParams;
	};

	FitThread(int _iterations, int _thread_id, int _n_threads, double** _data, double _length);
	~FitThread(void);

	boost::random::mt19937 * gen;
	boost::random::uniform_int_distribution<> * dist;

	UINT calculateMinT();
	static UINT startThread(LPVOID param);
	double randomDouble(double min, double max);
	double calculateResidual();

	threadContents contents;
	double** data;
	double length;
};

