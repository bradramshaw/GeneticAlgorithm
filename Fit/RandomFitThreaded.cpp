#include "StdAfx.h"
#include "RandomFitThreaded.h"
#include "ArrayAllocator.h"


#define cout std::cout
#define endl std::endl

#define PI 3.1415926535897932384626433832795028841971



RandomFitThreaded::RandomFitThreaded(double ** _data_set,int _length, int _n_threads, int _iterations)
{
	iterations = _iterations;
	min_parameters.A = 0;
	min_parameters.dF1 = 0;
	min_parameters.dF2 = 0;
	min_parameters.F = 0;
	min_parameters.m = 1.7;
	min_parameters.phi = 0;
	min_parameters.T = 4.2;
	min_parameters.Td = 0;
	min_parameters.chiSq = INFINITE;

	length = _length;
	data_set = _data_set;
	n_threads = _n_threads;

	p_threads = new FitThread*[n_threads];


}

RandomFitThreaded::~RandomFitThreaded(void)
{
	delete [] p_threads;
}

void RandomFitThreaded::setThreadNumber(int _n_threads)
{
	delete [] p_threads;
	n_threads = _n_threads;
	p_threads = new FitThread*[n_threads];
}

void RandomFitThreaded::setIterations(int _iterations)
{
	iterations = _iterations;
}



double RandomFitThreaded::getMinimum()
{
	double thread_iterations = (int) floor((double) iterations/n_threads);;

	for(int i = 0; i<n_threads; i++)
	{
		if(i == (n_threads - 1))
		{
			thread_iterations += (int) iterations%n_threads;
		}
		p_threads[i] = new FitThread(thread_iterations,i,n_threads,data_set,length);
	
	}


	HANDLE * handles = new HANDLE[n_threads];

	UINT n_threads_temp = n_threads;
	for(int i = 0; i < n_threads_temp; i++)
		{
			handles[i] = p_threads[i]->contents.hEvent;
		
			AfxBeginThread(p_threads[i]->startThread, reinterpret_cast<LPVOID>(p_threads[i]));			
		}
		
	WaitForMultipleObjects(n_threads_temp,handles,TRUE,INFINITE);
	
	for(int i = 0; i < n_threads_temp; i++)
	{
		if(p_threads[i]->contents.myParams.chiSq < min_parameters.chiSq)
		{
			min_parameters.A = p_threads[i]->contents.myParams.A;
			min_parameters.chiSq = p_threads[i]->contents.myParams.chiSq;
			min_parameters.dF1 = p_threads[i]->contents.myParams.dF1;
			min_parameters.dF2 = p_threads[i]->contents.myParams.dF2;
			min_parameters.F = p_threads[i]->contents.myParams.F;
			min_parameters.phi = p_threads[i]->contents.myParams.phi;
			min_parameters.Td = p_threads[i]->contents.myParams.Td;
		}
	}

	return 0;
}


void RandomFitThreaded::outputParams()
{
	cout<<"A: "<<min_parameters.A<<" "<<"F: "<<min_parameters.F<<" "<<"dF1: "<<min_parameters.dF1<<" "<<"dF2: "<<min_parameters.dF2<<" "<<"phi: "<<min_parameters.phi<<" "<<"Td: "<<min_parameters.Td<<" "<<endl<<"m: "<<min_parameters.m<<" "<<"T: "<<min_parameters.T<<endl<<"Residual: "<<min_parameters.chiSq<<endl<<endl;
	std::ofstream out;
	out.open("output.dat");
	out.precision(15);
	out<<min_parameters.A<<'\t'<<min_parameters.F<<'\t'<<min_parameters.dF1<<'\t'<<min_parameters.dF2<<'\t'<<min_parameters.phi<<'\t'<<min_parameters.Td<<'\t'<<min_parameters.m<<'\t'<<min_parameters.T;
	out.close();
}





