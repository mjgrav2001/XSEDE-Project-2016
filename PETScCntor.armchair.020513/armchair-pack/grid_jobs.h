//InitGridOffDiag
//InitGridDiag
//get_atom_col_block
//get_atom_connect_block
#ifndef GRID_JOBS_H_
#define GRID_JOBS_H_
#include <iostream>
#include <complex>
#include <new>
#include <ctime>
#include <fstream>
using namespace std;

//Returns the value of the block that an atom belongs to.
//For instance, in a 3600 atom model, atom 1359 will belong to block 133.
int get_atom_col_block(const PetscInt idim, int atom)
{
	return ((int) floor(atom / idim));
} //end get_atom_col_block


int get_atom_connect_block(const PetscInt idim, const PetscInt ndiag, int atom1, int atom2)
{
	int temp  = (int) floor(atom1 / idim);
	int temp2 = (int) floor(atom2 / idim);
	int result;

	if(temp == temp2)
	{
		//Then they are in the main diagonal
		result = 1;
	} //endif
	else if(temp == temp2+1 && temp2!= ndiag-1)
	{
		//They are in T/T-Dagger
		result = 0;
	} //endif
	else if(temp2 == temp+1 && temp != ndiag-1)
	{
		//They are in T/T-Dagger
		result = 2;
	} //endif
	else if(temp2 == 0 && temp == ndiag-1)
	{
		//They are in U/U-Dagger
		result = 3;
	} //endif
	else if(temp2 == ndiag-1 && temp == 0)
	{
		//They are in U/U-Dagger
		result = 4;
	} //endif
	else result = 999999;

	return (result);
} //end get_atom_connect_block

#endif /* GRID_JOBS_H_ */
