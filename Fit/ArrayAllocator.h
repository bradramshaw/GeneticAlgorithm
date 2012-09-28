#pragma once

namespace ArrayAllocator
{
	template<typename T>  T**  AllocateDynamicArray(int nRows, int nCols) //allocates an array
	{
		 T **dynamicArray;

		  dynamicArray = new T*[nRows];
		  for( int i = 0 ; i < nRows ; i++ )
		  dynamicArray[i] = new T [nCols];

		  return dynamicArray;
	} //allocates an array
	
	 template <typename T>  void  FreeDynamicArray(T** dArray,int nRows) // deallocates the same
	{
		for(int i = 0; i < nRows; i++)
			delete dArray[i];
		delete dArray;
	}
};
