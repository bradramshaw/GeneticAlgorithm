#include "stdafx.h"
#include "DataExtractor.h"
#include "ArrayAllocator.cpp"

double** dataArray;
int number_of_lines;

DataExtractor::DataExtractor(std::string name) // constructor
{
		std::vector<std::string> buffer1; //string buffer
		number_of_lines = dataReader(name, &buffer1);  //data reader puts the "name"d file into the buffer and returns its lenth in lines
		if(number_of_lines == 0)
			
			throw "File Not Found";  //empty file

		std::string ftype = ".dat";  //should be a .dat file
		char* p1 = (char*)name.c_str();  //converts name of file to a c string
		while(*p1!=NULL)
		{
			p1++; //goes to end of string
		}
		p1 = p1-4; //goes back 4 spaces
		if(ftype.compare(p1)!=0)
			throw "Wrong File Type"; //checks if the last 4 characters are .dat
		
		
		dataArray = ArrayAllocator::AllocateDynamicArray<double>(number_of_lines,4); //allocates an array of the appropriate size, right now 3 wide MAY NEED ADJUSTMENT
		bool y = dataConverter(&buffer1,dataArray);		
}

DataExtractor::~DataExtractor() // destructor 
{
	ArrayAllocator::FreeDynamicArray<double>(dataArray,number_of_lines);
}

int DataExtractor::dataReader(std::string name, std::vector<std::string>* buffer)	// Reads the data file into a vector in string format.
{	
	int x = 0; // return value
	std::string a;	//buffer
	
	std::ifstream in; // file to stream.
	in.open(name); //opens file
	if(in.is_open()) //checks for opened file
	{		
		while(getline(in,a))  //successivly gets lines from file (Each line has an x and a y data value)
		{
			(*buffer).push_back(a); //stores in the vector buffer
		}
		x=(*buffer).size();	// number of lines
	}
	in.close();
	return x; // returns if successful
}

bool DataExtractor::dataConverter(std::vector<std::string>* stringDat, double ** numDat )	// Converts a string format vector of numbers into a 2D array.
{		
	std::stringstream sbuff1,sbuff2; //buffers for string streams
	std::string buff1;	// stringn buffer	
	std::vector<std::string> vbuff1;	// vector buffer.

	for(int i=0; i < (int) (*stringDat).size(); i++) // for each line in the buffer (each data file line).
	{
		sbuff1 << (*stringDat)[i];	//put ith element into strinb buffer
		while(getline(sbuff1,buff1,'\t')) //get elements of each line (x and y)
		{
			vbuff1.push_back(buff1);	// put elements into buffer
		}		
		for(int j=0; j < 4; j++) //for x and y
		{ 
			double scale = 1;
			if(i == 0){
			scale = (1-0.0025);
			}
			numDat[i][j] =  scale * atof(vbuff1[j].c_str()); // convert string to float	and stores it
			
		} 
		vbuff1.clear(); //clears the buffers
		sbuff1.clear();
		sbuff2.clear();
		buff1.clear();
	}
	return true;
}

double** DataExtractor::getDataArray() //returns the data array
{
		return dataArray;
}

int DataExtractor::getNumberOfLines()
{
	return number_of_lines;
}

