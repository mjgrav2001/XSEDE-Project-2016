//InitGridOffDiag
//InitGridDiag
//get_atom_col_block
//get_atom_connect_block
#ifndef GRID_JOBS_H_
#define GRID_JOBS_H_
#include <iostream>
#include <complex>
using namespace std;

//Since the T,t-dagger,U, and U-dagger blocks are set once in the beginning, and then are never
//needed to be up dated, this only sets up those blocks and ignores the diagonal blocks.
void InitGridOffDiag(Parameters p, block stm[][3], double atomCoordinates[][3], double disCoordinates[][3])
{
	for(int i=0; i<p.ndiag-1; i++)
	{
		stm[i+1][0].set_values(p.rows, p.idim, p.ndiag, p.vhop, i, 0, p.energy, p.BField, 'T', atomCoordinates, disCoordinates);
	} //endfor

	for(int i=0; i<p.ndiag-1; i++)
	{
		stm[i][2].set_values(p.rows, p.idim, p.ndiag, p.vhop, i, 2, p.energy, p.BField, 'W', atomCoordinates, disCoordinates);
	} //endfor

	stm[0][0].set_values(p.rows, p.idim, p.ndiag, p.vhop, 0, 0, p.energy, p.BField, 'U', atomCoordinates, disCoordinates);

	stm[p.ndiag-1][2].set_values(p.rows, p.idim, p.ndiag, p.vhop, p.ndiag-1, 2, p.energy, p.BField, 'Q', atomCoordinates, disCoordinates);
} //end InitGridOffDiag


//Sets up only the diagonal blocks
void InitGridDiag(Parameters p, block stm[][3], double atomCoordinates[][3], double disCoordinates[][3])
{
	for(int i=0; i<p.ndiag; i++)
	{
		stm[i][1].set_values(p.rows, p.idim, p.ndiag, p.vhop, i, 1, p.energy, p.BField, 'D', atomCoordinates, disCoordinates);
	} //endfor
} //end InitGridDiag


//Returns the value of the block that an atom belongs to.
//For instance, in a 3600 atom model, atom 1359 will belong to block 133.
int get_atom_col_block(const PetscInt idim, int atom)
{
	return ((int) floor(atom / idim));
} //end get_atom_col_block


int get_atom_connect_block(const PetscInt idim, const PetscInt ndiag, int atom1, int atom2)
{
	int temp = (int) floor(atom1 / idim);
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
