// Fit.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "DataExtractor.h"
#include "GeneticAlgorithm2.h"

using namespace std;
int _tmain(int argc, _TCHAR* argv[])
{		
	LARGE_INTEGER time1,time2,freq;  // stores times and CPU frequency for profiling
	QueryPerformanceFrequency(&freq); //gets CPU frequency	
	double scaleFactor, crossingProbability; // how far the algorithm "looks" when mutating, and chance of mutation
	int nGenerations;  // number of generations, which is basically number of iterations
	try{
		DataExtractor extractor("C:/Users/root/Dropbox/Thesis/C++ data analysis/fullData.dat");  // get the data, which is field/amplitude/angle number. maybe bad format? should be consistent with mathematica maybe.
		double ** data  = extractor.getDataArray(); //puts data in 2x2 array of doubles.
		int nLines = extractor.getNumberOfLines(); // number of lines in the data set
		cout<<"Scale Factor: ";     
		cin>>scaleFactor;     // set scale factor. what is a good value? .1? i forget
		cout<<"Crossing Probability: ";
		cin>>crossingProbability;  // crossing prob. probably around .7 is good?
		GeneticAlgorithm2 geneticAlgorithm2(data,nLines,200,scaleFactor,crossingProbability);  // constructs the algorithm object. dataset/lines/population size/scale/crossing.

		while(true)      // population is set forever. maybe change this?
			{	
				cout<<"Number of generations: ";
				cin>>nGenerations;     // number of generationsn ie iterations. 
			QueryPerformanceCounter(&time1);	
				geneticAlgorithm2.calculateMinimum();
				geneticAlgorithm2.printMinimumParameters();	
				
				geneticAlgorithm2.calculateNewGenerations(nGenerations);		

				geneticAlgorithm2.calculateMinimum();
	
				geneticAlgorithm2.printMinimumParameters();		
			

			QueryPerformanceCounter(&time2);
			cout<<"Time for "<<"genetic algorithm"<<" is: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<endl<<endl;
			
			}	
		}
	catch(const char* str)
	{
		cout<<str<<endl;
	}
	return 0;
}

/*
{
	LARGE_INTEGER time1,time2,freq;
	QueryPerformanceFrequency(&freq);

	try{
		DataExtractor myExtractor("data2.dat");
		double ** myData  = myExtractor.getDataArray();
		int lines = myExtractor.getNumberOfLines();
		RandomFitThreadedFast fitter(myData,lines,4,100);
		int iterations,threads;
		
		while(true)
			{								
				cout<<"number of iterations: ";
				cin>>iterations;		
				cout<<"number of threads: ";
				cin>>threads;	
				cout<<endl;
			QueryPerformanceCounter(&time1);			
				fitter.setIterations(iterations);	
				fitter.setThreadNumber(threads);
				fitter.getMinimum();
				fitter.outputParams();	
			QueryPerformanceCounter(&time2);
			cout<<"Time for "<<threads<<" threads is: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<endl<<endl;

			}	
		}
	catch(const char* str)
	{
		cout<<str<<endl;
	}
	return 0;
}

/*
/*
int _tmain(int argc, _TCHAR* argv[])
{

	
	LARGE_INTEGER time1,time2,freq;
	QueryPerformanceFrequency(&freq);

	while(1)
	{
		int iterations,n_threads;
		cout<<"Iterations: ";
		cin>>iterations;
		cout<<endl;
		cout<<"Number of threads: ";
		cin>>n_threads;
		cout<<endl;

		double reg_total = 0;
		double t_total = 0;
							
		QueryPerformanceCounter(&time1);
		reg_total = sumNumbers(iterations);
		QueryPerformanceCounter(&time2);
		cout<<"unthreaded sum is: "<<reg_total<<endl;		
		cout<<"unthreaded time is: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<endl;
		cout<<endl;
						
		threadContents* threads = new threadContents[n_threads];
		HANDLE * handles = new HANDLE[n_threads];
		for(int i = 0; i < n_threads; i++)
		{
			threads[i].hEvent = CreateEvent(NULL, FALSE, FALSE, NULL);
			threads[i].total = 0;
			threads[i].iterations = iterations;
			threads[i].n_threads = n_threads;
			threads[i].thread_id = i;

			handles[i] = threads[i].hEvent;
		}
		

		QueryPerformanceCounter(&time1);

		for(int i = 0; i < n_threads; i++)
		{
			AfxBeginThread(sumNumbersT, &threads[i]);
		}
		
		WaitForMultipleObjects(n_threads,handles,TRUE,INFINITE);
		
		for(int i  = 0; i < n_threads; i++)
		{
			t_total += threads[i].total;
		}

		QueryPerformanceCounter(&time2);
		cout<<"threaded sum is: "<<t_total<<endl;		
		cout<<"threaded time is: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<endl<<endl;
		
		delete[] threads;
		delete[] handles;
	}

	return 0;
}

double sumNumbers(int _iterations)
{
	double total = 0;
	for(int i = 0; i<_iterations; i++)
	{
		total += sin(.2);
	}
	return total;
}

UINT sumNumbersT(LPVOID param)
{
	{
		
		threadContents* content_p = (threadContents*) param;
		int iterations = 0;
		if((*content_p).thread_id < ((*content_p).n_threads - 1))
		{
			iterations = (int) floor((double) (*content_p).iterations/(*content_p).n_threads);
		}
		else
		{
			iterations = (int) floor((double) (*content_p).iterations/(*content_p).n_threads)+ (int) (*content_p).iterations%(*content_p).n_threads;
		}
		

		for(int i = 0; i < iterations; i++)
		{
			(*content_p).total += sin(.2);
		}
		SetEvent((*content_p).hEvent);
		
	}

	return 0;
}
*/

//using namespace std;
//int _tmain(int argc, _TCHAR* argv[])
//{
//	try{
//		DataExtractor myExtractor("data2.dat");
//	
//		double ** myData  = myExtractor.getDataArray();
//		int lines = myExtractor.getNumberOfLines();
//		RandomFit myFitter(myData,lines);
//		int n;
//		
//		while(true)
//			{	
//							
//				cout<<"number of iterations: ";
//				cin>>n;				
//				
//				myFitter.calculateMin(n);
//				myFitter.outputParams();		
//			}	
//		ArrayAllocator::FreeDynamicArray<double>(myData,2);
//		}
//	catch(const char* str)
//	{
//		cout<<str<<endl;
//	}
//	return 0;
//}