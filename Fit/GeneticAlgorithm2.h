#pragma once
#include "Parameters2.h"

class GeneticAlgorithm2
{
public:
	GeneticAlgorithm2(double** dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability);
	~GeneticAlgorithm2(void);
	
	void calculateMinimum();
	void printMinimumParameters();
	void exportChiSq();
	void calculateNewGenerations(int nGenerations);
	
private:


	boost::random::mt19937 * _randomNumberGenerator;
	boost::random::uniform_int_distribution<> * _doubleDistribution;
	boost::random::uniform_int_distribution<> * _integerDistribution;
    VSLStreamStatePtr stream;
	int * ints1;
	int * ints2;
	int * ints3;

	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double ** _dataSet;
	
	double ** _residualArray;
	double ** _paramArray;
//  double ** _xVals;

	Parameters::fitParameters * _populationParametersOld, *_populationParametersNew;
	Parameters::fitParameters _minimumParameters;

	void initializeRandomNumberGenerators();
	void initializeParameters(double** dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability);
	double randomDouble(double min, double max);
	double calculateResidual(Parameters::fitParameters * parameters,int threadID);
	double calculateResidual2(Parameters::fitParameters * parameters,int threadID);
	double getYValue(double H, int theta,double T, double const * parameters);
	double getYValue2(int index, double H, int theta,double T, double const * parameters);
	double integrateLegendre( double H, int theta, double T, double * parameters);
	void resetParameters(int nPopulation, double scaleFactor, double crossingProbability);
	static	UINT startResidualThread(LPVOID param);
	void residualCalculatingThread(Parameters::arrayBounds * arrayBounds);

	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm2* pThis;
	};

};

