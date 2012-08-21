#include "StdAfx.h"
#include "FitThread.h"

#define PI 3.1415926535897932384626433832795028841971

#define cout std::cout
#define endl std::endl

FitThread::FitThread(int _iterations, int _thread_id, int _n_threads, double** _data, double _length)
{
	SYSTEMTIME t;
	GetLocalTime(&t);
	gen = new boost::random::mt19937(t.wMilliseconds+_thread_id);
	dist = new boost::random::uniform_int_distribution<>(1, RAND_MAX);

	contents.myParams.A = 0;
	contents.myParams.dF1 = 0;
	contents.myParams.dF2 = 0;
	contents.myParams.F = 0;
	contents.myParams.m = 1.7;
	contents.myParams.phi = 0;
	contents.myParams.T = 4.2;
	contents.myParams.Td = 0;
	contents.myParams.chiSq = INFINITE;

	contents.hEvent = CreateEvent(NULL, FALSE, FALSE, NULL);
	contents.myEnvel = new Envelope(0,0,0,0);
	contents.myInt = new Integrator();
	contents.myOscFun = new OscillationFunction(0,0,0,0);

	contents.n_threads = _n_threads;
	contents.thread_id = _thread_id;
	contents.iterations = _iterations;

	data = _data;
	length = _length;
}


FitThread::~FitThread(void)
{
}

UINT FitThread::startThread(LPVOID param)
{
	FitThread* _this = reinterpret_cast<FitThread *>(param);	
	return _this->calculateMinT();
}

UINT FitThread::calculateMinT()
{
	{		
			double F,dF1,dF2,phi,A,Td,T,m;	
			T = 4.2;
			m = 1.7;
			contents.myParams.m = m;
			contents.myParams.T = T;			

			for(int i = 0; i <contents.iterations; i++)
		{
			F = randomDouble(450,500);
			dF1 = randomDouble(30,50);
			dF2 = randomDouble(-20,0);
			phi = randomDouble(0,.5);
			A = randomDouble(0,2);
			Td = randomDouble(3,10);
	
			contents.myOscFun->changeParams(F,dF1,dF2,phi);
			contents.myEnvel->changeParams(A,T,Td,m);
				
			double residual = calculateResidual();
	
			if(residual < contents.myParams.chiSq)
			{
				contents.myParams.A = A;
				contents.myParams.chiSq = residual;
				contents.myParams.dF1 = dF1;
				contents.myParams.dF2 = dF2;
				contents.myParams.F = F;
				contents.myParams.phi = phi;
				contents.myParams.Td = Td;
			}	
		}	
	}
	SetEvent(contents.hEvent);
	return 0;}

double FitThread::randomDouble(double min, double max)
{				
	double randNum = min + (max - min)*( (double) (*dist)(*gen)/((double)RAND_MAX+1.0));
	return randNum;
}

double FitThread::calculateResidual()
{
	double residual = 0;

	for(int i  =0; i<length; i++)
	{		
		double env = contents.myEnvel->getYValue(data[i][0], 0);
		double osc = contents.myInt->integrateLegendre( contents.myOscFun,OscillationFunction::wrapper_getYValue, data[i][0],0,-PI,PI);
		residual += pow(env*osc - data[i][1],2);
	}		
	return residual;
}