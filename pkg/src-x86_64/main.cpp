// FloodFill-Algorithmus zum Finden der Hotspots
// Autor: Michael Bothmann
// Datum: Juli 2013

#include <stdio.h>
#include "FloodFill.h"
#include <vector>
#include <fstream>
#include <cstring>
#include <iostream>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace std;

int main(int argc, char *argv[])
{
	clock_t t1,t2;
	t1=clock();

	if (argc < 4)
	{
		cout << "Error, not enough parameters provided: $filename, $imageSize, $originalFarbe";
		return 0;
	}

	std::string fileName = argv[1];
	const int iImageBufferSize = 50000 + 2;

	const int iOriginalColor = 1;

	//std::ifstream file(fileName, ios::in|ios::binary);
	std::ifstream file(fileName.c_str(), ios::in|ios::binary);
	char* pImageBuffer = new char[iImageBufferSize * sizeof(double)];
	if (!file)
	{
		cout << "Error opening file\n";
	//	cout << strerror(errno);
		return 0;
	}


	file.read(pImageBuffer, iImageBufferSize * sizeof(double));
	
	file.close();

	int* pIntImageBuffer = reinterpret_cast<int*>(pImageBuffer);
	int iYDimension = reinterpret_cast<int*>(pImageBuffer)[0];
	int iXDimension = reinterpret_cast<int*>(pImageBuffer)[1];

	double* image = new double[iXDimension * iYDimension * 2];

	memcpy(image, &pIntImageBuffer[2], sizeof(double) * iXDimension * iYDimension);
	delete[] pImageBuffer;

	CFloodFill FloodFill(image, iXDimension, iYDimension);
	FloodFill.FloodCluster(iOriginalColor);

	//std::ofstream outputFile(fileName, ios::out|ios::binary);
	std::ofstream outputFile(fileName.c_str(), ios::out|ios::binary);
	outputFile.write(reinterpret_cast<char*>(&iXDimension), sizeof(int));
	outputFile.write(reinterpret_cast<char*>(&iYDimension), sizeof(int));
	outputFile.write(reinterpret_cast<char*>(image), sizeof(double) * 50000);
	outputFile.close();

	delete[] image;
	t2=clock();
    float diff ((float)t2-(float)t1);
	cout<< "Script runtime: " << diff << "ms" << endl;
	return 0;
}
