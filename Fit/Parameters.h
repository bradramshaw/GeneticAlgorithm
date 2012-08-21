#pragma once
namespace Parameters
{
	struct fitParameters
	{
		double A; //0
		double F; //1
		double dF1; //2
		double dF2; //3
		double Td; //4
		double phi; //5
		double ms; //6
		double chiSq; //7
		double T; //8
		double m; //9
	};

	struct arrayBounds
	{
		double start;
		double end;
		HANDLE handle;
		double time;
	};
}