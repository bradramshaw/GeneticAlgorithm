#include "StdAfx.h"
#include "GeneticAlgorithm.h"

static const double PI = 3.1415926535897932384626433832795028841971;
static const double TWOPI = 6.28318530717958647692528676656;
static const double OSC_CONST = 14.693185881921599;
static const double fToKfConstant = 0.0644735999999999; // converts a frequency to its kF value
static const double integrationNodes[] = {-0.999210123227436,-0.995840525118838,-0.989787895222222,-0.981067201752598,-0.969701788765053,-0.955722255839996,-0.939166276116423,-0.920078476177628,-0.898510310810046,-0.874519922646899,-0.84817198478593,-0.819537526162146,-0.788693739932264,-0.7557237753065849,-0.7207165133557299,-0.683766327381355,-0.644972828489477,-0.60444059704851,-0.562278900753944,-0.518601400058569,-0.473525841761707,-0.427173741583078,-0.379670056576798,-0.331142848268447,-0.281722937423261,-0.231543551376029,-0.180739964873425,-0.129449135396944,-0.07780933394953649,-0.0259597723012478,0.0259597723012475,0.07780933394953619,0.129449135396945,0.180739964873425,0.231543551376029,0.281722937423261,0.331142848268448,0.379670056576797,0.427173741583078,0.473525841761707,0.518601400058569,0.562278900753944,0.60444059704851,0.644972828489477,0.683766327381355,0.7207165133557299,0.7557237753065849,0.788693739932264,0.819537526162146,0.84817198478593,0.874519922646899,0.898510310810046,0.920078476177627,0.939166276116424,0.955722255839996,0.969701788765053,0.981067201752598,0.989787895222222,0.995840525118838,0.999210123227436};
static const double integrationWeights[] = {0.00202681196887378,0.00471272992695272,0.00738993116334546,0.0100475571822879,0.0126781664768157,0.015274618596785,0.0178299010142067,0.0203371207294572,0.0227895169439981,0.0251804776215212,0.027503556749925,0.0297524915007892,0.0319212190192966,0.03400389272494641,0.0359948980510843,0.0378888675692433,0.0396806954523818,0.04136555123558489,0.04293889283593489,0.0443964787957866,0.04573437971611471,0.0469489888489124,0.04803703181997171,0.04899557545575631,0.049822035690551,0.0505141845325093,0.0510701560698553,0.0514884515009809,0.05176794317491009,0.0519078776312198,0.0519078776312194,0.0517679431749094,0.0514884515009812,0.0510701560698546,0.05051418453250919,0.04982203569055041,0.0489955754557564,0.04803703181997189,0.0469489888489116,0.0457343797161143,0.04439647879578651,0.0429388928359353,0.041365551235585,0.0396806954523812,0.0378888675692443,0.0359948980510843,0.0340038927249456,0.031921219019296,0.0297524915007897,0.0275035567499253,0.0251804776215213,0.0227895169439967,0.0203371207294583,0.0178299010142076,0.0152746185967851,0.0126781664768161,0.0100475571822893,0.007389931163346141,0.004712729926953551,0.00202681196887367};
static const double besselCoeffs[] = {1.00000000001502,-0.0000000003971037387806575,-0.2500000025918557,0.00000004370805323228311,0.01562481926022373,0.0000003746964069017977,-0.0004344981963593569,0.0000003923846648059681,0.000006552041596216036,0.000000097577427433624,-0.0000000985582498286686,0.000000007259332734370126,-0.000000000815801230627031,0.0000000001695776376652846,-0.00000000001859116735850654,0.000000000001055434172075196,-0.00000000000003087777334068834,0.0000000000000003739549721015043};
static const int nBesselCoeffs = 18;
static const int nIntegrationNodes = 60;
static const int nIntegrationNodesP5 = 30;
static const double cosVals[] = {-0.999996921152254,-0.999914623060671,-0.999485408008636,-0.998231637389513,-0.995473361564615,-0.990340818503465,-0.981793088848555,-0.968644443742091,-0.949599907910554,-0.923301315943285,-0.888384632984825,-0.84354852950276,-0.78763315670609,-0.719706811490216,-0.63915679336933,-0.545779366724293,-0.439862510986826,-0.322254253103936,-0.194409017857362,-0.0584047664520222,0.0830751660609409,0.226799647599964,0.369088112552775,0.505947764490605,0.633244063417292,0.746895422959172,0.843080025699657,0.918440552064793,0.970271738184731,0.996676231418194,0.996676231418194,0.970271738184731,0.918440552064792,0.843080025699657,0.746895422959172,0.633244063417292,0.505947764490602,0.369088112552778,0.226799647599964,0.0830751660609409,-0.0584047664520222,-0.194409017857362,-0.322254253103936,-0.439862510986826,-0.545779366724293,-0.63915679336933,-0.719706811490216,-0.78763315670609,-0.84354852950276,-0.888384632984825,-0.923301315943285,-0.949599907910554,-0.96864444374209,-0.981793088848556,-0.990340818503465,-0.995473361564615,-0.998231637389513,-0.999485408008636,-0.999914623060671,-0.999996921152254};
//cosine of the kz values, which are the integration nodes (?)
static const double cos2Vals[] = {0.999987684627976,0.999658506821126,0.997942161644381,0.992932803770696,0.98193442716951,0.961549873588226,0.927835338621573,0.876544116784848,0.803479970207467,0.704970640045202,0.578454512247163,0.423148243252539,0.240731979085599,0.0359557890108273,-0.182957186979672,-0.404249765716061,-0.613041942856729,-0.792304392712849,-0.924410267551472,-0.993177766511369,-0.986197033567894,-0.897123839697065,-0.72754793034446,-0.488033719213919,-0.198003912293513,0.115705545674722,0.421567859467468,0.687066095354163,0.882854491840037,0.986727020547946,0.986727020547946,0.882854491840038,0.687066095354159,0.421567859467468,0.115705545674722,-0.198003912293513,-0.488033719213925,-0.727547930344456,-0.897123839697065,-0.986197033567894,-0.993177766511369,-0.924410267551472,-0.792304392712849,-0.613041942856729,-0.404249765716061,-0.182957186979672,0.0359557890108273,0.240731979085599,0.423148243252539,0.578454512247163,0.704970640045202,0.803479970207467,0.876544116784845,0.927835338621576,0.961549873588226,0.98193442716951,0.992932803770696,0.997942161644381,0.999658506821126,0.999987684627976};
//cosine of 2 times these values (second harmonic of warping)
static const int nVars = 7; // number of variables, F, dF, etc...
static const int nParams = 10;  //total parameters, including variables, but also T, m, and theta  which are set
static const double angles[] = {0.3869, 10.8832, 20.0916, -28.0054, 29.9469, 34.8886, 39.8032, 45.4099, 47.3163, 51.459, 54.4545, 57.3389};
//angles (Thesis version. offset of 1.06). not actually used, just the cos and tan of them. 
static const double tanAngles[] = {0.00675278133602108,0.192266475531181,0.365780915104845,-0.531830948943413,0.576115430708589,0.697313577268109,0.833264335024851,1.01441192849189,1.08430719141366,1.25532725297066,1.39959505375021,1.55999106515322};
// tangent of the angles (converted to radians of course)
static const double cosAngles[] = {0.999977200751846,0.982013952351089,0.939144893588774,0.88290311601117,0.866488343950678,0.820265824900911,0.76824734835429,0.702029882810786,0.677950975380776,0.623074545709823,0.581349533358446,0.539668305448932};
// cosine of the angles (converted to radians of course)

double totalTime = 0;  // for timing things
LARGE_INTEGER freq;


GeneticAlgorithm::GeneticAlgorithm(double** dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability){	
	QueryPerformanceFrequency(&freq);
	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	initializeParameters(dataSet, dataSetLength, nPopulation, scaleFactor, crossingProbability);	
}

GeneticAlgorithm::~GeneticAlgorithm(void){
	delete _doubleDistribution;
	delete _integerDistribution;
	delete _randomNumberGenerator;
	delete [] _populationParametersNew;
	delete [] _populationParametersOld;
}

void GeneticAlgorithm::initializeRandomNumberGenerators(){
	SYSTEMTIME t;
	GetLocalTime(&t);
	_randomNumberGenerator = new boost::random::mt19937(t.wMilliseconds);
	//_randomNumberGenerator = new boost::random::mt19937(0);
	_doubleDistribution = new boost::random::uniform_int_distribution<>(0, RAND_MAX);
	unsigned int max = _nPopulation - 1;
	_integerDistribution = new boost::random::uniform_int_distribution<>(0,max);
}


void GeneticAlgorithm::initializeParameters(double** dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability){

	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_dataSetLength = dataSetLength;
	_dataSet = dataSet;
	_populationParametersOld = new Parameters::fitParameters[nPopulation];
	_populationParametersNew = new Parameters::fitParameters[nPopulation];
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].A = randomDouble(0,5);
		_populationParametersOld[i].F = randomDouble(450,500);
		_populationParametersOld[i].dF1 = randomDouble(30,50);
		_populationParametersOld[i].dF2 = randomDouble(-20,0);
		_populationParametersOld[i].phi = randomDouble(0,.5);
		_populationParametersOld[i].Td = randomDouble(3,10);
		_populationParametersOld[i].ms = randomDouble(1.0,2.2);
		_populationParametersOld[i].T = 4.2;
		_populationParametersOld[i].m = 1.7;
		_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i]);
	}

	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersNew[i].T = 4.2;
		_populationParametersNew[i].m = 1.7;
	}
		_minimumParameters.A = 1;
		_minimumParameters.F = 0;
		_minimumParameters.dF1 = 0;
		_minimumParameters.dF2 = 0;
		_minimumParameters.phi = 0;
		_minimumParameters.Td = 0;
		_minimumParameters.ms = 0;
		_minimumParameters.T = 4.2;
		_minimumParameters.m = 1.7;
		_minimumParameters.chiSq = INFINITE;
}

void GeneticAlgorithm::resetParameters(int nPopulation, double scaleFactor, double crossingProbability){
	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_nPopulation = nPopulation;
	delete [] _populationParametersOld;
	delete [] _populationParametersNew;
	_populationParametersOld = new Parameters::fitParameters[nPopulation];
	_populationParametersNew = new Parameters::fitParameters[nPopulation];
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].A = randomDouble(0,5);
		_populationParametersOld[i].F = randomDouble(450,510);
		_populationParametersOld[i].dF1 = randomDouble(10,50);
		_populationParametersOld[i].dF2 = randomDouble(-20,0);
		_populationParametersOld[i].phi = randomDouble(0,.5);
		_populationParametersOld[i].Td = randomDouble(3,10);
		_populationParametersOld[i].ms = randomDouble(1.0,2.2);
		_populationParametersOld[i].T = 4.2;
		_populationParametersOld[i].m = 1.7;
		
	}
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i]);
	}

	delete _integerDistribution;
	_integerDistribution = new boost::random::uniform_int_distribution<>(0, _nPopulation-1);

}

double GeneticAlgorithm::randomDouble(double min, double max){				
	double randNum = min + (max - min)*( (double) (*_doubleDistribution)(*_randomNumberGenerator)/((double)RAND_MAX));
	return randNum;
}

void GeneticAlgorithm::calculateMinimum(){
	for(int i  = 0; i < _nPopulation; i++){
		if(_populationParametersOld[i].chiSq < _minimumParameters.chiSq){
			_minimumParameters.A = _populationParametersOld[i].A;
			_minimumParameters.F = _populationParametersOld[i].F;
			_minimumParameters.dF1 = _populationParametersOld[i].dF1;
			_minimumParameters.dF2 = _populationParametersOld[i].dF2;
			_minimumParameters.phi = _populationParametersOld[i].phi;
			_minimumParameters.Td = _populationParametersOld[i].Td;
			_minimumParameters.ms = _populationParametersOld[i].ms;
			_minimumParameters.chiSq = _populationParametersOld[i].chiSq;
		}
	}
}

void GeneticAlgorithm::printMinimumParameters(){
	std::cout<<"A: "<<_minimumParameters.A<<" "<<"F: "<<_minimumParameters.F<<" "<<"dF1: "<<_minimumParameters.dF1<<" "<<"dF2: "<<_minimumParameters.dF2<<" "<<"phi: "<<_minimumParameters.phi<<" "<<"Td: "<<_minimumParameters.Td<<" "<<std::endl<<"ms: "<<_minimumParameters.ms<<" "<<"m: "<<_minimumParameters.m<<" "<<"T: "<<_minimumParameters.T<<std::endl<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
		std::ofstream out;
		out.open("output.dat");
		out.precision(15);
		out<<_minimumParameters.A<<'\t'<<_minimumParameters.F<<'\t'<<_minimumParameters.dF1<<'\t'<<_minimumParameters.dF2<<'\t'<<_minimumParameters.phi<<'\t'<<_minimumParameters.Td<<'\t'<<_minimumParameters.ms<<'\t'<<_minimumParameters.m<<'\t'<<_minimumParameters.T;
		out.close();
}

double GeneticAlgorithm::calculateResidual(Parameters::fitParameters * parameters){
	/*double residual = 0;
	double modelDataDifference;
	#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){	
		modelDataDifference = integrateLegendre(_dataSet[i][0],_dataSet[i][2], parameters) - _dataSet[i][1];
		residual += modelDataDifference*modelDataDifference;
	}		
	return residual;*/

	double * residualArray = new double[_dataSetLength];
	double total = 0;
	double modelDataDifference;
	double * paramArray = new double[nParams];
	double * paramPointer = &(parameters->A);
	for(int i = 0; i < nParams; i++){
		paramArray[i] = *paramPointer;
		paramPointer++;
	}
#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){	
		modelDataDifference = integrateLegendre(_dataSet[i][0],_dataSet[i][2], paramArray) - _dataSet[i][1];
		residualArray[i] = modelDataDifference*modelDataDifference;
	}		
#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){
		total += residualArray[i];
	}
	delete [] residualArray;
	delete [] paramArray;
	return total;
}

double GeneticAlgorithm::integrateLegendre( double H, int theta, double * parameters){
	double total1 = 0;
	double kf_sqF_tan = fToKfConstant * sqrt(parameters[1]) * tanAngles[theta];
	double inv_H_cos = 1/(H*cosAngles[theta]);

#pragma ivdep
	for(int i = nIntegrationNodesP5; i < nIntegrationNodes; i++){		
	total1 += integrationWeights[i]*sin(TWOPI*((parameters[1] + (parameters[2])*cosVals[i]*boost::math::cyl_bessel_j(0,kf_sqF_tan)+ (parameters[3])*cos2Vals[i]*boost::math::cyl_bessel_j(0,2*kf_sqF_tan)  )*inv_H_cos-parameters[5]));
	//	total1 += integrationWeights[i]*getYValue(i,H,theta, parameters);
	}	
	//return sqrt(H)*(parameters[0])*TWOPI*(total1)* cos(TWOPI*parameters[6]/(2*cosAngles[theta]))*exp(-OSC_CONST*(parameters[4])*(parameters[9])*inv_H_cos)*(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)/sinh(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos);
	return cos(TWOPI*parameters[6]/(2*cosAngles[theta]))*exp(-OSC_CONST*(parameters[4])*(parameters[9])*inv_H_cos)*(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)/sinh(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)*sqrt(H)*(parameters[0])*TWOPI*(total1);
}

double GeneticAlgorithm::getYValue(int index, double H, int theta, double * parameters){	
	double kf_sqF_tan = fToKfConstant * sqrt(parameters[1]) * tanAngles[theta];
	double inv_H_cos = 1/(H*cosAngles[theta]);
	//return cos(TWOPI*parameters->ms/(2*cosAngles[theta]))*exp(-OSC_CONST*parameters->Td*parameters->m*inv_H_cos)*(OSC_CONST*parameters->T*parameters->m*inv_H_cos)/sinh(OSC_CONST*parameters->T*parameters->m*inv_H_cos)*sin(TWOPI*((parameters->F+parameters->dF1*cosVals[index]*besselJ(kf_sqF_tan) + parameters->dF2*cos2Vals[index]*besselJ(2*kf_sqF_tan))*inv_H_cos-parameters->phi));
	return cos(TWOPI*parameters[6]/(2*cosAngles[theta]))*exp(-OSC_CONST*(parameters[4])*(parameters[9])*inv_H_cos)*(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)/sinh(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)*sin(TWOPI*((parameters[1] + (parameters[3])*cos2Vals[index]*boost::math::cyl_bessel_j(0,2*kf_sqF_tan) + (parameters[2])*cosVals[index]*boost::math::cyl_bessel_j(0,kf_sqF_tan) )*inv_H_cos-parameters[5]));
	//double arg1 = (parameters[3])*cos2Vals[index]*boost::math::cyl_bessel_j(0,2*kf_sqF_tan);
	//double arg2 = (parameters[2])*cosVals[index]*boost::math::cyl_bessel_j(0,kf_sqF_tan);
	//double arg3 = TWOPI*((parameters[1] + arg1 + arg2 )*inv_H_cos-parameters[5]);
	//double result1 = sin(arg3);
	//double result2 = cos(TWOPI*parameters[6]/(2*cosAngles[theta]))*exp(-OSC_CONST*(parameters[4])*(parameters[9])*inv_H_cos);
	//double result3 = (OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos)/sinh(OSC_CONST*(parameters[8])*(parameters[9])*inv_H_cos);
	//return result1*result2*result3;
}

 double GeneticAlgorithm::besselJ(double x){
	double total = 0;
	double xPow = 1;
	int i = 0;	
	do
	{
		total += besselCoeffs[i]*xPow;
		xPow *= x;
		i++;
	}
	while(i < nBesselCoeffs);

	return total;
}

void GeneticAlgorithm::calculateNewGenerations(int nGenerations){
	double averageTime = 0;

	for(int i = 0; i < nGenerations; i++){
	
		for(int j = 0; j < _nPopulation; j++){
			int g1 = (*_integerDistribution)(*_randomNumberGenerator);
			int g2 = (*_integerDistribution)(*_randomNumberGenerator);
			int g3 = (*_integerDistribution)(*_randomNumberGenerator);
			double * pointerToOldVariable = &_populationParametersOld[j].A;
			double * pointerToNewVariable = &_populationParametersNew[j].A;
			double * pointerTog1Variable = &_populationParametersOld[g1].A;
			double * pointerTog2Variable = &_populationParametersOld[g2].A;
			double * pointerTog3Variable = &_populationParametersOld[g3].A;
		
			for(int k = 0; k < nVars; k++){
				double p = randomDouble(0,1);
				if(p > _crossingProbability){
					*pointerToNewVariable = *pointerToOldVariable;
				}
				else{
					*pointerToNewVariable = *pointerTog1Variable + _scaleFactor*(*pointerTog2Variable - *pointerTog3Variable);
				}
				pointerToNewVariable++;
				pointerToOldVariable++;
				pointerTog1Variable++;
				pointerTog2Variable++;
				pointerTog3Variable++;
			}		
		}	
		totalTime = 0;

		const UINT nThreads = 4;
		HANDLE threadEvents[nThreads];
		Parameters::arrayBounds threadBounds[nThreads];
		threadContents threadContents[nThreads];
		for(int m = 0; m<nThreads; m++){
			threadEvents[m] = CreateEvent(NULL, FALSE, FALSE, NULL);
			int nPopulationPerThread =  _nPopulation/nThreads;
			if(m != (nThreads-1)){
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1;
			}
			else{
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1 + _nPopulation%nThreads;
			}
			threadBounds[m].handle = threadEvents[m];
			threadBounds[m].time = 0;
			threadContents[m].arrayBounds = threadBounds[m];
			threadContents[m].pThis = this;
			AfxBeginThread(startResidualThread, (LPVOID) &threadContents[m]);		
			
		}
		WaitForMultipleObjects(nThreads,threadEvents,TRUE,INFINITE);			
		
		totalTime = threadContents[0].arrayBounds.time+threadContents[1].arrayBounds.time+threadContents[2].arrayBounds.time+threadContents[3].arrayBounds.time;
		averageTime +=totalTime;
		
		/*for(int l = 0; l < _nPopulation; l++){
			_populationParametersNew[l].chiSq = calculateResidual(&_populationParametersNew[l]);
			if(_populationParametersNew[l].chiSq < _populationParametersOld[l].chiSq){
				_populationParametersOld[l] = _populationParametersNew[l];
			}
		}*/
		/*calculateMinimum();
		std::ofstream out;
		out.open("output_c_p1.dat",std::ios::app);
		out.precision(15);
		out<<_minimumParameters.chiSq<<std::endl;
		out.close();*/
	}
	std::cout<<"Average time per generation per thread is: "<<averageTime/(4*nGenerations)<<"ms"<<std::endl<<std::endl;
}

UINT GeneticAlgorithm::startResidualThread(LPVOID param){
	threadContents * contents = (threadContents*) param;
	contents->pThis->residualCalculatingThread(&(contents->arrayBounds));
	return 0;
}

void GeneticAlgorithm::residualCalculatingThread(Parameters::arrayBounds * arrayBounds){
	LARGE_INTEGER time1,time2;
	#pragma ivdep
	for(int i = arrayBounds->start; i <= arrayBounds->end; i++){
		QueryPerformanceCounter(&time1);
			_populationParametersNew[i].chiSq = calculateResidual(&_populationParametersNew[i]);
		QueryPerformanceCounter(&time2);
		arrayBounds->time += 1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart);
			if(_populationParametersNew[i].chiSq < _populationParametersOld[i].chiSq){
				_populationParametersOld[i] = _populationParametersNew[i];
			}
	}
	SetEvent(arrayBounds->handle);
}