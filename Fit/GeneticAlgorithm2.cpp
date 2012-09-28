#include "StdAfx.h"
#include "GeneticAlgorithm2.h"

static const double PI = 3.1415926535897932384626433832795028841971;
static const double TWOPI = 6.28318530717958647692528676656;
static const double OSC_CONST = 14.693185881921599;
static const double fToKfConstant = 0.0647245399554491; // converts a frequency to its kF value
static const double integrationNodes[] = {-0.999210123227436,-0.995840525118838,-0.989787895222222,-0.981067201752598,-0.969701788765053,-0.955722255839996,-0.939166276116423,-0.920078476177628,-0.898510310810046,-0.874519922646899,-0.84817198478593,-0.819537526162146,-0.788693739932264,-0.7557237753065849,-0.7207165133557299,-0.683766327381355,-0.644972828489477,-0.60444059704851,-0.562278900753944,-0.518601400058569,-0.473525841761707,-0.427173741583078,-0.379670056576798,-0.331142848268447,-0.281722937423261,-0.231543551376029,-0.180739964873425,-0.129449135396944,-0.07780933394953649,-0.0259597723012478,0.0259597723012475,0.07780933394953619,0.129449135396945,0.180739964873425,0.231543551376029,0.281722937423261,0.331142848268448,0.379670056576797,0.427173741583078,0.473525841761707,0.518601400058569,0.562278900753944,0.60444059704851,0.644972828489477,0.683766327381355,0.7207165133557299,0.7557237753065849,0.788693739932264,0.819537526162146,0.84817198478593,0.874519922646899,0.898510310810046,0.920078476177627,0.939166276116424,0.955722255839996,0.969701788765053,0.981067201752598,0.989787895222222,0.995840525118838,0.999210123227436};
static const double integrationWeights[] = {0.00202681196887378,0.00471272992695272,0.00738993116334546,0.0100475571822879,0.0126781664768157,0.015274618596785,0.0178299010142067,0.0203371207294572,0.0227895169439981,0.0251804776215212,0.027503556749925,0.0297524915007892,0.0319212190192966,0.03400389272494641,0.0359948980510843,0.0378888675692433,0.0396806954523818,0.04136555123558489,0.04293889283593489,0.0443964787957866,0.04573437971611471,0.0469489888489124,0.04803703181997171,0.04899557545575631,0.049822035690551,0.0505141845325093,0.0510701560698553,0.0514884515009809,0.05176794317491009,0.0519078776312198,0.0519078776312194,0.0517679431749094,0.0514884515009812,0.0510701560698546,0.05051418453250919,0.04982203569055041,0.0489955754557564,0.04803703181997189,0.0469489888489116,0.0457343797161143,0.04439647879578651,0.0429388928359353,0.041365551235585,0.0396806954523812,0.0378888675692443,0.0359948980510843,0.0340038927249456,0.031921219019296,0.0297524915007897,0.0275035567499253,0.0251804776215213,0.0227895169439967,0.0203371207294583,0.0178299010142076,0.0152746185967851,0.0126781664768161,0.0100475571822893,0.007389931163346141,0.004712729926953551,0.00202681196887367};
static const int nIntegrationNodes = 60;
static const int nIntegrationNodesP5 = 30;

static const double cosVals[] = {-0.999996921152254,-0.999914623060671,-0.999485408008636,-0.998231637389513,-0.995473361564615,-0.990340818503465,-0.981793088848555,-0.968644443742091,-0.949599907910554,-0.923301315943285,-0.888384632984825,-0.84354852950276,-0.78763315670609,-0.719706811490216,-0.63915679336933,-0.545779366724293,-0.439862510986826,-0.322254253103936,-0.194409017857362,-0.0584047664520222,0.0830751660609409,0.226799647599964,0.369088112552775,0.505947764490605,0.633244063417292,0.746895422959172,0.843080025699657,0.918440552064793,0.970271738184731,0.996676231418194,0.996676231418194,0.970271738184731,0.918440552064792,0.843080025699657,0.746895422959172,0.633244063417292,0.505947764490602,0.369088112552778,0.226799647599964,0.0830751660609409,-0.0584047664520222,-0.194409017857362,-0.322254253103936,-0.439862510986826,-0.545779366724293,-0.63915679336933,-0.719706811490216,-0.78763315670609,-0.84354852950276,-0.888384632984825,-0.923301315943285,-0.949599907910554,-0.96864444374209,-0.981793088848556,-0.990340818503465,-0.995473361564615,-0.998231637389513,-0.999485408008636,-0.999914623060671,-0.999996921152254};
//cosine of the kz values, which are the integration nodes (?)
static const double cos2Vals[] = {0.999987684627976,0.999658506821126,0.997942161644381,0.992932803770696,0.98193442716951,0.961549873588226,0.927835338621573,0.876544116784848,0.803479970207467,0.704970640045202,0.578454512247163,0.423148243252539,0.240731979085599,0.0359557890108273,-0.182957186979672,-0.404249765716061,-0.613041942856729,-0.792304392712849,-0.924410267551472,-0.993177766511369,-0.986197033567894,-0.897123839697065,-0.72754793034446,-0.488033719213919,-0.198003912293513,0.115705545674722,0.421567859467468,0.687066095354163,0.882854491840037,0.986727020547946,0.986727020547946,0.882854491840038,0.687066095354159,0.421567859467468,0.115705545674722,-0.198003912293513,-0.488033719213925,-0.727547930344456,-0.897123839697065,-0.986197033567894,-0.993177766511369,-0.924410267551472,-0.792304392712849,-0.613041942856729,-0.404249765716061,-0.182957186979672,0.0359557890108273,0.240731979085599,0.423148243252539,0.578454512247163,0.704970640045202,0.803479970207467,0.876544116784845,0.927835338621576,0.961549873588226,0.98193442716951,0.992932803770696,0.997942161644381,0.999658506821126,0.999987684627976};
//cosine of 2 times these values (second harmonic of warping)
static const int nVars = 14; // number of variables, F, dF, etc...
static const int nParams = 16;  //total parameters, including variables, but also T, and theta  which are set

static const double angles[] = {0.00675268,0.189949,0.350664, 0.488787,0.522672,0.608921,0.694698,0.792552,0.825825,0.898129,0.95041,1.00075};
//angles (Thesis version. offset of 1.06). not actually used, just the cos and tan of them. 
static const double tanAngles[] = {0.00675278133602108,0.192266475531181,0.365780915104845,-0.531830948943413,0.576115430708589,0.697313577268109,0.833264335024851,1.01441192849189,1.08430719141366,1.25532725297066,1.39959505375021,1.55999106515322};
// tangent of the angles (converted to radians of course)
static const double cosAngles[] = {0.999977200751846,0.982013952351089,0.939144893588774,0.88290311601117,0.866488343950678,0.820265824900911,0.76824734835429,0.702029882810786,0.677950975380776,0.623074545709823,0.581349533358446,0.539668305448932};
// cosine of the angles (converted to radians of course)
static const double sinAngles[] = {0.00675262737768367,0.188808361540989,0.343521278592944,-0.469555202013317,0.499197305479117,0.571982496672431,0.640153115861043,0.712147487281022,0.73510711803128,0.782162457861855,0.813653931388473,0.841877734646712};
static const int nThreads = 4;
double totalTime = 0;  // for timing things
LARGE_INTEGER freq;
	_CrtMemState s1,s2,s3;
/* upper bounds
static const double angles[] = {0.4089,10.9056,20.1194,-27.9829,29.9799,34.9083,39.8237,45.4472,47.3352,51.479,54.362,57.3587};
static const double tanAngles[] = {0.00713677,0.192671,0.366331,-0.531326,0.576882,0.697824,0.833868,1.01573,1.08503,1.25623,1.39483,1.56118};
static const double cosAngles[] = {0.999975,0.98194,0.938978,0.883088,0.866201,0.820069,0.768019,0.701566,0.677708,0.622802,0.582663,0.539377};
static const double sinAngles[] = {0.00713659,0.189191,0.343977,-0.469208,0.499696,0.572264,0.640427,0.712605,0.735331,0.78238,0.812714,0.842064};*/

// lower bounds
//static const double angles[] = {0.3649,10.8616,20.0754,-28.0269,29.9359,34.8643,39.7797,45.4032,47.2912,51.435,54.318,57.3147};
//static const double tanAngles[] = {0.00636879,0.191874,0.365461,-0.532312,0.575859,0.696683,0.832567,1.01418,1.08335,1.25425,1.39257,1.55854};
//static const double cosAngles[] = {0.99998,0.982085,0.939242,0.882727,0.866584,0.820509,0.768511,0.702113,0.678273,0.623402,0.583287,0.540024};
//static const double sinAngles[] = {0.00636866,0.188437,0.343256,-0.469886,0.499031,0.571634,0.639837,0.712066,0.73481,0.781901,0.812266,0.84165};

GeneticAlgorithm2::GeneticAlgorithm2(double** dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability){	
	QueryPerformanceFrequency(&freq);
	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	initializeParameters(dataSet, dataSetLength, nPopulation, scaleFactor, crossingProbability);	
}

GeneticAlgorithm2::~GeneticAlgorithm2(void){
//	delete _doubleDistribution;
	//delete _integerDistribution;
//	delete _randomNumberGenerator;
	delete [] _populationParametersNew;
	delete [] _populationParametersOld;
}

void GeneticAlgorithm2::initializeRandomNumberGenerators(){
	SYSTEMTIME t;
	GetLocalTime(&t);
	//_randomNumberGenerator = new boost::random::mt19937(t.wMilliseconds);
	//_randomNumberGenerator = new boost::random::mt19937(0);
//	_doubleDistribution = new boost::random::uniform_int_distribution<>(0, RAND_MAX);
	unsigned int max = _nPopulation - 1;
	// _integerDistribution = new boost::random::uniform_int_distribution<>(0,max);

	vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
}


void GeneticAlgorithm2::initializeParameters(double** dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability){

	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_dataSetLength = dataSetLength;
	_dataSet = dataSet;

	_residualArray = new double*[nThreads];
	_paramArray = new double*[nThreads];
	for(int i = 0; i < nThreads; i++){
	_residualArray[i] = new double[_dataSetLength];
	_paramArray[i] = new double[nParams];
	// xVals[i] = new double[4];
	}

	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].A1 = randomDouble(0.0,1.0);
		_populationParametersOld[i].A2 = randomDouble(0.0,1.0);
		_populationParametersOld[i].F1 = randomDouble(460,480);
		_populationParametersOld[i].F2 = randomDouble(520,540);
		_populationParametersOld[i].dF1 = randomDouble(25.,45.0);
		_populationParametersOld[i].dF2 = 0;
		_populationParametersOld[i].phi1 = randomDouble(0.0,1.0);
		_populationParametersOld[i].phi2 = randomDouble(0.0,1.0);
		_populationParametersOld[i].Td1 = randomDouble(4., 8.);
		_populationParametersOld[i].Td2 = randomDouble(4., 8.);
	//	_populationParametersOld[i].ms1 = randomDouble(1,2.2);
	//	_populationParametersOld[i].ms2 = randomDouble(1.5,1.7);
		_populationParametersOld[i].ms1 = randomDouble(1.0,1.4);
		_populationParametersOld[i].ms2 = randomDouble(1.5,1.7);
	//	_populationParametersOld[i].m1 =  randomDouble(1.0,3);
	//	_populationParametersOld[i].m2 =  randomDouble(1.0,3);
		_populationParametersOld[i].m1 = 1.7;
		_populationParametersOld[i].m2 =  1.7;
	//	_populationParametersOld[i].dF12 = randomDouble(-15,5);
	//	_populationParametersOld[i].ms11 = randomDouble(0.0,.2);
	//	_populationParametersOld[i].ms22 = randomDouble(0.0,.2);
		_populationParametersOld[i].T = 4.2;
		_populationParametersOld[i].chiSq = calculateResidual2(&_populationParametersOld[i],0);
	}

	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersNew[i].T = 0;
		_populationParametersNew[i].m1 = 0;
		_populationParametersNew[i].m2 = 0;
	}
		_minimumParameters.A1 = 1;
		_minimumParameters.A2 = 1;
		_minimumParameters.F1 = 0;
		_minimumParameters.F2 = 0;
		_minimumParameters.dF1 = 0;
		_minimumParameters.dF2 = 0;
		_minimumParameters.phi1 = 0;
		_minimumParameters.phi2 = 0;
		_minimumParameters.Td1 = 0;
		_minimumParameters.Td2 = 0;
		_minimumParameters.ms1 = 0;
		_minimumParameters.ms2 = 0;
		_minimumParameters.m1 = 0;
		_minimumParameters.m2 = 0;
	//	_minimumParameters.dF12 = 0;
	//	_minimumParameters.ms11 = 0;
	//	_minimumParameters.ms22 = 0;
		_minimumParameters.T = 0;
		_minimumParameters.chiSq = INFINITE;
}

void GeneticAlgorithm2::resetParameters(int nPopulation, double scaleFactor, double crossingProbability){
	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_nPopulation = nPopulation;
	delete [] _populationParametersOld;
	delete [] _populationParametersNew;
//	_populationParametersOld = new Parameters::fitParameters[nPopulation];
//	_populationParametersNew = new Parameters::fitParameters[nPopulation];
	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].A1 = randomDouble(0.0,5.0);
		_populationParametersOld[i].A2 = randomDouble(0.0,5.0);
		_populationParametersOld[i].F1 = randomDouble(440.0,500.0);
		_populationParametersOld[i].F2 = randomDouble(500.0,560.0);
		_populationParametersOld[i].dF1 = randomDouble(20.0,60.0);
		_populationParametersOld[i].dF2 =0;
		_populationParametersOld[i].phi1 = randomDouble(0.0,1.0);
		_populationParametersOld[i].phi2 = randomDouble(0.0,1.0);
		_populationParametersOld[i].Td1 = randomDouble(2.0,10.0);
		_populationParametersOld[i].Td2 = randomDouble(2.0,10.0);
		_populationParametersOld[i].ms1 = randomDouble(1,2.6);
		_populationParametersOld[i].ms2 = randomDouble(1,2.6);
	  //  _populationParametersOld[i].m1 =  randomDouble(1.3,1.9);
	//	_populationParametersOld[i].m2 =  randomDouble(1.3,1.9);
		_populationParametersOld[i].m1 = 1.7;
		_populationParametersOld[i].m2 =  1.7;
		_populationParametersOld[i].dF12 =  randomDouble(-15,5);
	//	_populationParametersOld[i].ms11 = randomDouble(0.0,.2);
	//	_populationParametersOld[i].ms22 = randomDouble(0.0,.2);
		_populationParametersOld[i].T = 4.2;
		
	}
	for(int i  = 0; i < _nPopulation; i++){
		_populationParametersOld[i].chiSq = calculateResidual2(&_populationParametersOld[i],0);
	}

	//delete _integerDistribution;
	delete [] ints1;
	delete [] ints2;
	delete [] ints3;
//	_integerDistribution = new boost::random::uniform_int_distribution<>(0, _nPopulation-1);
	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
	
}

double GeneticAlgorithm2::randomDouble(double min, double max){				
	double randNum;
	 vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &randNum, min, max);
	return randNum;
}

void GeneticAlgorithm2::calculateMinimum(){
	for(int i  = 0; i < _nPopulation; i++){
		if(_populationParametersOld[i].chiSq < _minimumParameters.chiSq){
			_minimumParameters.A1 = _populationParametersOld[i].A1;
			_minimumParameters.A2 = _populationParametersOld[i].A2;
			_minimumParameters.F1 = _populationParametersOld[i].F1;
			_minimumParameters.F2 = _populationParametersOld[i].F2;
			_minimumParameters.dF1 = _populationParametersOld[i].dF1;
			_minimumParameters.dF2 = _populationParametersOld[i].dF2;
			_minimumParameters.phi1 = _populationParametersOld[i].phi1;
			_minimumParameters.phi2 = _populationParametersOld[i].phi2;
			_minimumParameters.Td1 = _populationParametersOld[i].Td1;
			_minimumParameters.Td2 = _populationParametersOld[i].Td2;
			_minimumParameters.ms1 = _populationParametersOld[i].ms1;
			_minimumParameters.ms2 = _populationParametersOld[i].ms2;
			_minimumParameters.m1 = _populationParametersOld[i].m1;
			_minimumParameters.m2 = _populationParametersOld[i].m2;
		//	_minimumParameters.dF12 = _populationParametersOld[i].dF12;
		//	_minimumParameters.ms11 = _populationParametersOld[i].ms11;
		//	_minimumParameters.ms22 = _populationParametersOld[i].ms22;
			_minimumParameters.chiSq = _populationParametersOld[i].chiSq;
		}
	}
}

void GeneticAlgorithm2::printMinimumParameters(){
	std::cout<<"A1: "<<_minimumParameters.A1<<" "<<"A2: "<<_minimumParameters.A2<<" "<<"F1: "<<_minimumParameters.F1<<" "<<"F2: "<<_minimumParameters.F2<<" "<<"dF1: "<<_minimumParameters.dF1<<" "<<"dF2: "<<_minimumParameters.dF2<<" "<<std::endl<<"phi1: "<<_minimumParameters.phi1<<" "<<"phi2: "<<_minimumParameters.phi2<<" "<<"Td1: "<<_minimumParameters.Td1<<" "<<"Td2: "<<_minimumParameters.Td2<<" "<<std::endl<<"ms1: "<<_minimumParameters.ms1<<" "<<"ms2: "<<_minimumParameters.ms2<<" "<<"m1: "<<_minimumParameters.m1<<" "<<"m2: "<<_minimumParameters.m2<<std::endl<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
		std::ofstream out;
		out.open("output.dat");
		out.precision(15);
		out<<_minimumParameters.A1<<'\t'<<_minimumParameters.A2<<'\t'<<_minimumParameters.F1<<'\t'<<_minimumParameters.F2<<'\t'<<_minimumParameters.dF1<<'\t'<<_minimumParameters.dF2<<'\t'<<_minimumParameters.phi1<<'\t'<<_minimumParameters.phi2<<'\t'<<_minimumParameters.Td1<<'\t'<<_minimumParameters.Td2<<'\t'<<_minimumParameters.ms1<<'\t'<<_minimumParameters.ms2<<'\t'<<_minimumParameters.m1<<'\t'<<_minimumParameters.m2;
		out.close();
}

double GeneticAlgorithm2::calculateResidual(Parameters::fitParameters * parameters, int threadID){
	
	double total = 0;	
	double * paramPointer = &(parameters->A1);

	for(int i = 0; i < nParams; i++){
		_paramArray[threadID][i] = *paramPointer;
		paramPointer++;
	}

#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){

		_residualArray[threadID][i] = getYValue(_dataSet[i][0],_dataSet[i][2], _dataSet[i][3],_paramArray[threadID]) - _dataSet[i][1];
	}	

#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){
		total += _residualArray[threadID][i]*_residualArray[threadID][i];
	}

	return total;
}

double GeneticAlgorithm2::calculateResidual2(Parameters::fitParameters * parameters, int threadID){
	
	double total = 0;	
	double * paramPointer = &(parameters->A1);

	for(int i = 0; i < nParams; i++){
		_paramArray[threadID][i] = *paramPointer;
		paramPointer++;
	}

#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){

		_residualArray[threadID][i] = integrateLegendre(_dataSet[i][0],_dataSet[i][2], _dataSet[i][3],_paramArray[threadID]) - _dataSet[i][1];
	}	

#pragma ivdep
	for(int i = 0; i < _dataSetLength; i++){
		total += _residualArray[threadID][i]*_residualArray[threadID][i];
	}

	return total;
}

double GeneticAlgorithm2::integrateLegendre( double H, int theta, double T, double * parameters){
	double total1 = 0;
	double total2 = 0;
	double kf_sqF1_tan = fToKfConstant * sqrt(parameters[2]) * tanAngles[theta];
	double kf_sqF1_tan2 = 2.0*fToKfConstant * sqrt(parameters[2]) * tanAngles[theta];
	double kf_sqF2_tan = fToKfConstant * sqrt(parameters[3]) * tanAngles[theta];
	double inv_H_cos = 1.0/(H*cosAngles[theta]);;

#pragma ivdep
	for(int i = nIntegrationNodesP5; i < nIntegrationNodes; i++){		
	//	+ (parameters[14])*cos2Vals[i]*j0(kf_sqF1_tan2) ) + atan((1+OSC_CONST*(parameters[6])*(parameters[12])*inv_H_cos)/(PI *(parameters[4])*j0(kf_sqF1_tan)*inv_H_cos ))
	total1 += TWOPI * integrationWeights[i]*cos(TWOPI*((parameters[2] + (parameters[4])*cosVals[i]*j0(kf_sqF1_tan) )*inv_H_cos+parameters[8]));
	total2 += TWOPI * integrationWeights[i]*cos(TWOPI*((parameters[3] + (parameters[5])*cosVals[i]*j0(kf_sqF2_tan) )*inv_H_cos+parameters[9]));
	}	
	//return parameters[0]*cos(TWOPI*(parameters[10]*sqrt(cosAngles[theta]*cosAngles[theta] + (1.0 - parameters[15])*sinAngles[theta]*sinAngles[theta] ))/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[6])*(parameters[12])*inv_H_cos)*(OSC_CONST*T*(parameters[12])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[12])*inv_H_cos)*total1 + parameters[1]*cos(TWOPI*(parameters[11]*sqrt(cosAngles[theta]*cosAngles[theta] + (1.0 - parameters[16])*sinAngles[theta]*sinAngles[theta] ))/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[7])*(parameters[13])*inv_H_cos)*(OSC_CONST*T*(parameters[13])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[13])*inv_H_cos)*total2;
	return parameters[0]*cos(TWOPI*(parameters[10])/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[6])*(parameters[12])*inv_H_cos)*(OSC_CONST*T*(parameters[12])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[12])*inv_H_cos)*total1 + parameters[1]*cos(TWOPI*(parameters[11])/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[7])*(parameters[13])*inv_H_cos)*(OSC_CONST*T*(parameters[13])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[13])*inv_H_cos)*total2;

}

double GeneticAlgorithm2::getYValue( double H, int theta, double T, double const * parameters){	

	double kf_sqF1_tan = fToKfConstant * sqrt(parameters[2]) * tanAngles[theta];
	double kf_sqF2_tan = fToKfConstant * sqrt(parameters[3]) * tanAngles[theta];
	double inv_H_cos = 1.0/(H*cosAngles[theta]);

double surf1 = parameters[0]*cos(TWOPI*(parameters[10]+(H/30.) * parameters[15])/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[6])*(parameters[12])*inv_H_cos)*(OSC_CONST*T*(parameters[12])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[12])*inv_H_cos)*cos(TWOPI*(parameters[2]*inv_H_cos+parameters[8]))*j0(TWOPI* parameters[4]*inv_H_cos *j0( kf_sqF1_tan ) );
double surf2 = parameters[1]*cos(TWOPI*(parameters[11]+(H/30.) * parameters[16])/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[7])*(parameters[13])*inv_H_cos)*(OSC_CONST*T*(parameters[13])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[13])*inv_H_cos)*cos(TWOPI*(parameters[3]*inv_H_cos+parameters[9]))*j0(TWOPI* parameters[5]*inv_H_cos *j0( kf_sqF2_tan ) );

//double surf1 = parameters[0]*cos(TWOPI*parameters[10]/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[6])*(parameters[12])*inv_H_cos * pow(1.0 + 0.296461556242335 * H*H * cosAngles[theta]*cosAngles[theta] / ((parameters[6])*(parameters[12])*(parameters[6])*(parameters[12])),.25))*(OSC_CONST*T*(parameters[12])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[12])*inv_H_cos)*cos(TWOPI*(parameters[2]*inv_H_cos+parameters[8]))*j0(TWOPI* parameters[4]*inv_H_cos *j0( kf_sqF1_tan ) );
//double surf2 = parameters[1]*cos(TWOPI*parameters[11]/(2.0*cosAngles[theta]))*exp(-OSC_CONST*(parameters[7])*(parameters[13])*inv_H_cos * pow(1.0 + 0.296461556242335 * H*H * cosAngles[theta]*cosAngles[theta] / ((parameters[7])*(parameters[13])*(parameters[7])*(parameters[13])),.25))*(OSC_CONST*T*(parameters[13])*inv_H_cos)/sinh(OSC_CONST*T*(parameters[13])*inv_H_cos)*cos(TWOPI*(parameters[3]*inv_H_cos+parameters[9]))*j0(TWOPI* parameters[5]*inv_H_cos *j0( kf_sqF2_tan ) );

return (surf1+surf2);
}


void GeneticAlgorithm2::calculateNewGenerations(int nGenerations){
	double averageTime = 0;

	for(int i = 0; i < nGenerations; i++){


	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints1, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints2, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints3, 0, _nPopulation );
	 				
		for(int j = 0; j < _nPopulation; j++){
		//	int g1 = (*_integerDistribution)(*_randomNumberGenerator);
		//	int g2 = (*_integerDistribution)(*_randomNumberGenerator);
		//	int g3 = (*_integerDistribution)(*_randomNumberGenerator);
		
			double * pointerToOldVariable = &_populationParametersOld[j].A1;
			double * pointerToNewVariable = &_populationParametersNew[j].A1;
			double * pointerTog1Variable = &_populationParametersOld[ints1[j]].A1;
			double * pointerTog2Variable = &_populationParametersOld[ints2[j]].A1;
			double * pointerTog3Variable = &_populationParametersOld[ints3[j]].A1;
		
	//could be optimzed for vector arithmetic
			for(int k = 0; k < nVars; k++){
				double p = randomDouble(0,1);
				if(p > _crossingProbability){
					*pointerToNewVariable = *pointerToOldVariable;
				}
				else{
					*pointerToNewVariable = *pointerTog1Variable + _scaleFactor*(*pointerTog2Variable - *pointerTog3Variable);
					if((k == 8) || (k == 9))
						if(*pointerToNewVariable > 1.0)
							*pointerToNewVariable = *pointerToNewVariable - 1.0;
						if(*pointerToNewVariable < -1.0)
							*pointerToNewVariable = *pointerToNewVariable + 1.0;
				}
				pointerToNewVariable++;
				pointerToOldVariable++;
				pointerTog1Variable++;
				pointerTog2Variable++;
				pointerTog3Variable++;
			}		
		}	

		totalTime = 0;

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
			threadBounds[m].threadID = m;
			threadContents[m].arrayBounds = threadBounds[m];
			threadContents[m].pThis = this;

			AfxBeginThread(startResidualThread, (LPVOID) &threadContents[m]);		
		}
	
		//_CrtMemCheckpoint( &s2 );		
//		_CrtMemDumpStatistics( &s2 );
	//	std::cout<<s2.lTotalCount<<std::endl;

		WaitForMultipleObjects(nThreads,threadEvents,TRUE,INFINITE);	


		totalTime = threadContents[0].arrayBounds.time+threadContents[1].arrayBounds.time+threadContents[2].arrayBounds.time+threadContents[3].arrayBounds.time;
		averageTime +=totalTime;
	//	calculateMinimum();
	//	exportChiSq();
	}
	std::cout<<"Average time per generation per thread is: "<<averageTime/(4*nGenerations)<<"ms"<<std::endl<<std::endl;
	
}

UINT GeneticAlgorithm2::startResidualThread(LPVOID param){
	threadContents * contents = (threadContents*) param;
	contents->pThis->residualCalculatingThread(&(contents->arrayBounds));
	return 0;
}

void GeneticAlgorithm2::residualCalculatingThread(Parameters::arrayBounds * arrayBounds){
	LARGE_INTEGER time1,time2;
	#pragma ivdep
	for(int i = arrayBounds->start; i <= arrayBounds->end; i++){
		QueryPerformanceCounter(&time1);
		_populationParametersNew[i].chiSq = calculateResidual2(&_populationParametersNew[i],arrayBounds->threadID);
		QueryPerformanceCounter(&time2);
		arrayBounds->time += 1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart);
			if(_populationParametersNew[i].chiSq < _populationParametersOld[i].chiSq){
				_populationParametersOld[i] = _populationParametersNew[i];
			}
	}
	SetEvent(arrayBounds->handle);
	return;
}

void GeneticAlgorithm2::exportChiSq(){
    	std::ofstream out;
		out.open("chiSq10.dat",std::ios_base::app);
		out.precision(15);
		out<<_minimumParameters.chiSq<<"\t";
		out.close();
}