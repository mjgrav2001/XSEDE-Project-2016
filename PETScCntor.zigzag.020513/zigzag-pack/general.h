//Contains general functions
//ReadInCoordinates
//ReadInDisplaced
//permuteRow
//permuteColumn
#ifndef GENERAL_H_
#define GENERAL_H_
#include <iostream>
#include <complex>
#include <new>
#include <ctime>
#include <fstream>
using namespace std;

void ReadInCoordinates(const PetscInt rows, double atomCoordinates[][3])
{
	ifstream inFile("ARRAY_posr.dat");
	//ifstream inFile("positions");
	int atomNum;

	if(inFile.is_open())
	{
		for(int lineNum=1; (lineNum<=rows && !inFile.eof()); lineNum++)
		{
			inFile >> atomCoordinates[lineNum-1][0] >> atomCoordinates[lineNum-1][1] >> atomCoordinates[lineNum-1][2] >> atomNum;
		} //endfor
		inFile.clear();              // forget we hit the end of file
		inFile.seekg(0, ios::beg); 
		inFile.close();
	} //endif
	else cout << "Error: Positions file not opened!\n";
} //end ReadInCoordinates


void ReadInDisplaced(const PetscInt rows, double disCoordinates[][3])
{
	ifstream inFile("ARRAY_posr_new.dat");
	//ifstream inFile("positions");
	int disNum;
	
	if(inFile.is_open())
	{
		for(int lineNum=1; (lineNum<=rows && !inFile.eof()); lineNum++)
		{
			inFile >> disCoordinates[lineNum-1][0] >> disCoordinates[lineNum-1][1] >> disCoordinates[lineNum-1][2] >> disNum;
		} //endfor
		inFile.clear();              // forget we hit the end of file
		inFile.seekg(0, ios::beg); 
		inFile.close();
	} //endif
	else cout << "Error: Displaced positions file not opened!\n";
} //end ReadInDisplaced


//******************************************************************************************************
void permuteRow (int leadsize, complex<double> **array, int rowA, int rowB)
{
	complex<double> tempVar;
	for(int j=0; j<leadsize; j++)
	{
		tempVar = array[rowA][j];
		array[rowA][j] = array[rowB][j];
		array[rowB][j] = tempVar;
	}
} //end permuteRow


//******************************************************************************************************
void permuteColumn (int leadsize, complex<double> **array, int colA, int colB)
{
	complex<double> tempVar;
	for(int i=0; i<leadsize; i++)
	{
		tempVar = array[i][colA];
		array[i][colA] = array[i][colB];
		array[i][colB] = tempVar;
	}
} //end permuteColumn

#endif /* GENERAL_H_ */
