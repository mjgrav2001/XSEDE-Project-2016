//Contains general functions
//ReadInCoordinates
//ReadInDisplaced
#ifndef GENERAL_H_
#define GENERAL_H_
#include <iostream>
#include <fstream>
#include <complex>
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

#endif /* GENERAL_H_ */
