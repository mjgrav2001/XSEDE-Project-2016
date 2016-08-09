//InitGridOffDiag
//InitGridDiag
//get_atom_col_block
//get_atom_connect_block
#ifndef GRID_JOBS_H_
#define GRID_JOBS_H_
#include <iostream>
#include <complex>
using namespace std;

//Returns the value of the block that an atom belongs to.
//For instance, in a 3600 atom model, atom 1359 will belong to block 133.
int get_atom_col_block(const PetscInt ndiag, int atom)
{
	return ((int) floor(atom / ndiag));
} //end get_atom_col_block


int get_atom_connect_block(const PetscInt idim, const PetscInt ndiag, int atom1, int atom2)
{
	int temp  = (int) floor(atom1 / ndiag);
	int temp2 = (int) floor(atom2 / ndiag);
	int result;

	if(temp == temp2)
	{
		//Then they are in the main diagonal
		result = 1;
	} //endif
	else if(temp == temp2+1 && temp2!= idim-1)
	{
		//They are in T/T-Dagger
		result = 0;
	} //endif
	else if(temp2 == temp+1 && temp != idim-1)
	{
		//They are in T/T-Dagger
		result = 2;
	} //endif
	else if(temp2 == 0 && temp == idim-1)
	{
		//They are in U/U-Dagger
		result = 3;
	} //endif
	else if(temp2 == idim-1 && temp == 0)
	{
		//They are in U/U-Dagger
		result = 4;
	} //endif
	else result = 999999;

	return (result);
} //end get_atom_connect_block

#endif /* GRID_JOBS_H_ */
