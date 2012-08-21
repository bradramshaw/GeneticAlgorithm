#pragma once
#include "Parameters.h"

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(double** dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability);
	~GeneticAlgorithm(void);
	
	void calculateMinimum();
	void printMinimumParameters();
	void calculateNewGenerations(int nGenerations);
	
private:

	boost::random::mt19937 * _randomNumberGenerator;
	boost::random::uniform_int_distribution<> * _doubleDistribution;
	boost::random::uniform_int_distribution<> * _integerDistribution;

	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double ** _dataSet;
	
	Parameters::fitParameters * _populationParametersOld, *_populationParametersNew;
	Parameters::fitParameters _minimumParameters;

	void initializeRandomNumberGenerators();
	void initializeParameters(double** dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability);
	double randomDouble(double min, double max);
	double calculateResidual(Parameters::fitParameters * parameters);
	double integrateLegendre( double H, int theta, double * parameters);
	double getYValue(int index, double H, int theta, double * parameters);
	double besselJ(double x);
	void resetParameters(int nPopulation, double scaleFactor, double crossingProbability);
	static	UINT startResidualThread(LPVOID param);
	void residualCalculatingThread(Parameters::arrayBounds * arrayBounds);

	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm* pThis;
	};

};

