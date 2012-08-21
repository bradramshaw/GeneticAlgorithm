#pragma once
namespace Parameters
{
	__declspec( align( 128) ) struct fitParameters
	{
		double A1; //0
		double A2;  //1
		double F1; //2
		double F2; //3
		double dF1; //4
		double dF2; //5
		double Td1; //6
		double Td2;  //7
		double phi1; //8
		double phi2; //9
		double ms1; //10
		double ms2; //11
		double m1; //12
		double m2; //13
		double dF12; //14
		double chiSq; //15
		double T; //16
	};

    struct arrayBounds
	{
		double start;
		double end;
		HANDLE handle;
		double time;
		int threadID;
	};
}