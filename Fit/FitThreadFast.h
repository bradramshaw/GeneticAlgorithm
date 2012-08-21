#pragma once
#include "Envelope.h"
#include "OscillationFunction.h"
#include "IntegratorFast.h"
#include "Parameters.h"

class FitThreadFast
{
public:
	struct threadContents
	{
		HANDLE hEvent;
		int n_threads;
		int thread_id;
		int iterations;
		IntegratorFast *myInt;
		Parameters::fitParameters myParams;
	};

	FitThreadFast(int _iterations, int _thread_id, int _n_threads, double** _data, double _length);
	~FitThreadFast(void);

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

