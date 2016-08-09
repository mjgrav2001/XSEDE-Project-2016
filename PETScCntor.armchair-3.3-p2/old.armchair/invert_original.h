#ifndef INVERT_H_
#define INVERT_H_
#include "petscksp.h"
#include "classes.h"
#include "general.h"
#include "grid_jobs.h"
#include <complex>
#include <iostream>
#include <ctime>
#include "petscksp.h"

using namespace std;

const double PI = 3.1415926535;
void permuteRow(complex<double> [4][4], int, int);
void permuteColumn(complex<double> [4][4], int, int);

#undef __FUNCT__
#define __FUNCT__ "InvertHamiltonian"
PetscErrorCode InvertHamiltonian (MPI_Comm PETSC_COMM_SUBGROUP, PetscInt size, PetscInt rank, PetscInt sizesg, PetscInt ranksg, Parameters& p, Mat hamiltonian, Mat inverseHamiltonian, Mat identityMat, PC pc, KSP ksp, 
								  int isl180[], int isr180[], int isl[], int isr[], complex<double> gl[4][4], complex<double> gr[4][4], 
								  double atomCoordinates[][3], double disCoordinates[][3], block blockHamiltonian)
{
	//Declarations
	time_t programStart;
	time(&programStart);
	
	//PETSc Setup
	PetscErrorCode ierr;
	
	//Variable Declarations';
	time(&p.invertStart);
	bool check;
	PetscInt idim1 = p.idim;
	PetscInt idim2 = idim1*idim1;
	PetscInt itor = p.idim/2;
	PetscInt itor2 = itor*itor;
	
    int ici,ici1;
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
	
	//PetscInt *globalRowIDs, *globalColIDs, *globalRowID, *globalColID;	
	//ierr = PetscMalloc(idim1*sizeof(PetscInt), &globalRowIDs);CHKERRQ(ierr);
	//ierr = PetscMalloc(idim1*sizeof(PetscInt), &globalColIDs);CHKERRQ(ierr);	
	//ierr = PetscMalloc(itor*sizeof(PetscInt), &globalRowID);CHKERRQ(ierr);
	//ierr = PetscMalloc(itor*sizeof(PetscInt), &globalColID);CHKERRQ(ierr);
	
	//PetscScalar *ham, *ham0;	
	//ierr = PetscMalloc(itor*itor*sizeof(PetscScalar), &ham0);CHKERRQ(ierr);		
	//ierr = PetscMalloc(idim1*idim1*sizeof(PetscScalar), &ham);CHKERRQ(ierr);
	
	int m = (int) p.ndiag / p.size0;

	//******************************************************************************************************
	for(int j=ranksg*m; j<(ranksg+1)*m; j++)
	{
		for(int j=0; j<p.ndiag; j++)
		{
			blockHamiltonian.set_zeros(p.idim, p.ndiag);
			blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j, 1, p.energy, p.BField, 'D', atomCoordinates, disCoordinates);
			blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
										 j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);			
			for(int k=0; k<p.idim; k++)
			{
				for(int l=0; l<p.idim; l++)
				{
					ham = blockHamiltonian.get_data(k,l);
					if(ham != complex<double>(0,0))
						ierr = MatSetValue(hamiltonian, k+p.idim*j, l+p.idim*j, ham, INSERT_VALUES);CHKERRQ(ierr);
				} //endfor
			} //endfor
			//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
			//							j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);				
		} //endfor
	}
	
	if (p.ndiag % p.size0 != 0)
	{
		if (ranksg==0)
		{
			for(int j=m*p.size0; j<p.ndiag; j++)
			{
				for(int j=0; j<p.ndiag; j++)
				{
					blockHamiltonian.set_zeros(p.idim, p.ndiag);
					blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j, 1, p.energy, p.BField, 'D', atomCoordinates, disCoordinates);
					blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
												 j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);			
					for(int k=0; k<p.idim; k++)
					{
						for(int l=0; l<p.idim; l++)
						{
							ham = blockHamiltonian.get_data(k,l);
							if(ham != complex<double>(0,0))
								ierr = MatSetValue(hamiltonian, k+p.idim*j, l+p.idim*j, ham, INSERT_VALUES);CHKERRQ(ierr);
						} //endfor
					} //endfor
					//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
					//							j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
					
				} //endfor
			}
		}
	}
	
	//Print run time
	if(rank == 0)
	{
		int programTime = (int) difftime(time(NULL), programStart);
		int seconds = programTime % 60;
		int minutes = (programTime / 60) % 60;
		int hours = programTime / 3600;
		cout << "\nTotal run time: ";
		if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
		else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
		else cout << seconds << " seconds" << endl;
	}
	
	//******************************************************************************************************
	for(int j=ranksg*m; j<(ranksg+1)*m; j++)
	{
		for(int j=0; j<p.ndiag-1; j++)
		{
			blockHamiltonian.set_zeros(p.idim, p.ndiag);
			blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j, 2, p.energy, p.BField, 'W', atomCoordinates, disCoordinates);
			blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
										 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
			for(int k=0; k<p.idim; k++)
			{
				for(int l=0; l<p.idim; l++)
				{
					ham = blockHamiltonian.get_data(k,l);
					if(ham != complex<double>(0,0))
						ierr = MatSetValue(hamiltonian, k+p.idim*j, l+p.idim*(j-1), ham, INSERT_VALUES);CHKERRQ(ierr);
					
					if(j == p.ndiag-2)
					{
						check = true;
						for(int i=1; i<=itor; i++)
						{
							if(check == true)
							{
								ici = 2*(i-1);
								ici1 = ici+1;
								ham = blockHamiltonian.get_data(ici,ici1);
								if(ham != complex<double>(0,0))
									ierr = MatSetValue(hamiltonian, p.rows-p.idim+ici, p.rows-2*p.idim+ici1, ham, INSERT_VALUES);CHKERRQ(ierr);
								check = !check;
							} //endif
							else
							{
								ici = 2*i-1;
								ici1 = ici-1;
								ham = blockHamiltonian.get_data(ici,ici1);
								if(ham != complex<double>(0,0))
									ierr = MatSetValue(hamiltonian, p.rows-p.idim+ici, p.rows-2*p.idim+ici1, ham, INSERT_VALUES);CHKERRQ(ierr);
								check = !check;
							} //endelse
						} //endfor
					} //endif
				} //endfor
			} //endfor
			//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
			//							 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);				
		} //endfor
	}
	
	if (p.ndiag % p.size0 != 0)
	{
		if (ranksg==0)
		{
			for(int j=m*p.size0; j<p.ndiag; j++)
			{
				for(int j=0; j<p.ndiag-1; j++)
				{
					blockHamiltonian.set_zeros(p.idim, p.ndiag);
					blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j, 2, p.energy, p.BField, 'W', atomCoordinates, disCoordinates);
					blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
												 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
					for(int k=0; k<p.idim; k++)
					{
						for(int l=0; l<p.idim; l++)
						{
							ham = blockHamiltonian.get_data(k,l);
							if(ham != complex<double>(0,0))
								ierr = MatSetValue(hamiltonian, k+p.idim*j, l+p.idim*(j-1), ham, INSERT_VALUES);CHKERRQ(ierr);
							
							if(j == p.ndiag-2)
							{
								check = true;
								for(int i=1; i<=itor; i++)
								{
									if(check == true)
									{
										ici = 2*(i-1);
										ici1 = ici+1;
										ham = blockHamiltonian.get_data(ici,ici1);
										if(ham != complex<double>(0,0))
											ierr = MatSetValue(hamiltonian, p.rows-p.idim+ici, p.rows-2*p.idim+ici1, ham, INSERT_VALUES);CHKERRQ(ierr);
										check = !check;
									} //endif
									else
									{
										ici = 2*i-1;
										ici1 = ici-1;
										ham = blockHamiltonian.get_data(ici,ici1);
										if(ham != complex<double>(0,0))
											ierr = MatSetValue(hamiltonian, p.rows-p.idim+ici, p.rows-2*p.idim+ici1, ham, INSERT_VALUES);CHKERRQ(ierr);
										check = !check;
									} //endelse
								} //endfor
							} //endif
						} //endfor
					} //endfor
					//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
					//							 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
				} //endfor
			}
		}
	}

	//******************************************************************************************************
	for(int j=ranksg*m; j<(ranksg+1)*m; j++)
	{
		blockHamiltonian.set_zeros(p.idim, p.ndiag);
		blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j-1, 0, p.energy, p.BField, 'T', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
									 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);	
		for(int k=0; k<p.idim; k++)
		{
			for(int l=0; l<p.idim; l++)
			{
				ham = blockHamiltonian.get_data(k,l);
				if(ham != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k+p.idim*(j-1), l+p.idim*j, ham, INSERT_VALUES);CHKERRQ(ierr);
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);			
	}
	
	if (p.ndiag % p.size0 != 0)
	{
		if (ranksg==0)
		{
			for(int j=m*p.size0; j<p.ndiag; j++)
			{
				blockHamiltonian.set_zeros(p.idim, p.ndiag);
				blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, j-1, 0, p.energy, p.BField, 'T', atomCoordinates, disCoordinates);
				blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
											 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);	
				for(int k=0; k<p.idim; k++)
				{
					for(int l=0; l<p.idim; l++)
					{
						ham = blockHamiltonian.get_data(k,l);
						if(ham != complex<double>(0,0))
							ierr = MatSetValue(hamiltonian, k+p.idim*(j-1), l+p.idim*j, ham, INSERT_VALUES);CHKERRQ(ierr);
					} //endfor
				} //endfor
				//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
				//							 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);		
			}
		}
	}
	
	//******************************************************************************************************
	if (ranksg==0)
	{
		blockHamiltonian.set_zeros(p.idim, p.ndiag);
		blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, 0, 0, p.energy, p.BField, 'U', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
									 0, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);		
		for(int k=0; k<p.idim; k++)
		{
			for (int l=0; l<p.idim; l++)
			{
				ham = blockHamiltonian.get_data(k,l);
				if(ham != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k, l+p.idim*(p.ndiag-1), ham, INSERT_VALUES);CHKERRQ(ierr);			
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 0, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	

		blockHamiltonian.set_zeros(p.idim, p.ndiag);
		blockHamiltonian.set_values(p.rows, p.idim, p.ndiag, p.vhop, p.ndiag-1, 2, p.energy, p.BField, 'Q', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
									 p.ndiag-1, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);		
		for(int k=0; k<p.idim; k++)
		{
			for (int l=0; l<p.idim; l++)
			{
				ham = blockHamiltonian.get_data(k,l);
				if(ham != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k+p.idim*(p.ndiag-1), l, ham, INSERT_VALUES);CHKERRQ(ierr);
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.idim, p.ndiag, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 p.ndiag-1, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
	}

	//******************************************************************************************************
	//Print run time
	if(rank == 0)
	{
		int programTime = (int) difftime(time(NULL), programStart);
		int seconds = programTime % 60;
		int minutes = (programTime / 60) % 60;
		int hours = programTime / 3600;
		cout << "\nTotal run time: ";
		if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
		else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
		else cout << seconds << " seconds" << endl;
	}
	
	//******************************************************************************************************
	//Assemble PETSc Hamiltonian
	ierr = MatAssemblyBegin(hamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(hamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//Reset PETSc Solver
	ierr = KSPSetOperators(ksp, hamiltonian, hamiltonian, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);
	
	//Factor and solve hamiltonian * inverseHamiltonian = identityMat for inverseHamiltonian
	ierr = PCFactorGetMatrix(pc, &hamiltonian);CHKERRQ(ierr);
	ierr = MatMatSolve(hamiltonian, identityMat, inverseHamiltonian);CHKERRQ(ierr);
	//******************************************************************************************************
	
	//Print run time
	if(rank == 0)
	{
		int programTime = (int) difftime(time(NULL), programStart);
		int seconds = programTime % 60;
		int minutes = (programTime / 60) % 60;
		int hours = programTime / 3600;
		cout << "\nTotal run time: ";
		if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
		else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
		else cout << seconds << " seconds" << endl;
	}
	
	//ierr = PetscFree(ham);CHKERRQ(ierr);
	//ierr = PetscFree(ham0);CHKERRQ(ierr);
	//ierr = PetscFree(globalRowIDs);CHKERRQ(ierr);
	//ierr = PetscFree(globalColIDs);CHKERRQ(ierr);
	//ierr = PetscFree(globalRowID);CHKERRQ(ierr);
	//ierr = PetscFree(globalColID);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
} //end InvertHamiltonian


#undef __FUNCT__
#define __FUNCT__ "RetrieveInformation"
PetscErrorCode RetrieveInformation (MPI_Comm PETSC_COMM_SUBGROUP, PetscInt size, PetscInt rank, PetscInt sizesg, PetscInt ranksg, Parameters p, Mat inverseHamiltonian, double IVvalues[], FILE* outFile, int isl[], int isr[], complex<double> gl[4][4], complex<double> gr[4][4])
{
	//Declarations
	time_t programStart;
	time(&programStart);

	//PETSc Setup
	const PetscScalar *rowValues;
	
	PetscScalar trace;
	PetscInt nonZeros;
	PetscErrorCode ierr;
	
	PetscInt numSubMatrices = 2;
	
	static Mat *localSubMatrices;
	static IS *indexSetRows, *indexSetColumns;
	
	int n,size0n;
	n =  (int) size / p.size0;
	size0n = (int) p.size0 * n;
	
	//Variable Declarations
	complex<double> gaml[4][4], gamrr[4][4];
	complex<double> prodl[4][4], prodr[4][4], zres[4][4];
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
	double densityOfStates = 0.0, therm = .03;
	double emu1 = -p.vsdMax/2.0, emu2 = p.vsdMax/2.0;
	int vsl[4],vsr[4];
	int currentStep = (p.energy < 0.0 ? (p.energy - p.energyMin) / p.energyStep + .5 : (p.energy - p.energyMin) / p.energyStep - .5);
	
	//cout <<  " rank " << rank << " size " << size << endl;
	//cout <<  " ranksg " << ranksg << " sizesg " << sizesg << endl;
	
	ierr = MatGetTrace(inverseHamiltonian, &trace);CHKERRQ(ierr);
	densityOfStates = -1.0 * trace.imag() / PI;

	if(p.size0 > 1) //parallel
	{
		//This only works if the rows and columns of gaml and gamrr are determined by arrays
		//Will not work for arbitrary gaml[i][j] = inverseHamiltonian[m][n]
		//if(energy == energyMin)
		//{
		
		for(int i=0; i<p.leadSize; i++)
		{
			vsl[i]=isl[i]-1;
			vsr[i]=isr[i]-1;
			for(int j=0; j<p.leadSize; j++)
			{
				gaml[i][j] = complex<double>(0.,0.);
				gamrr[i][j] = complex<double>(0.,0.);
			}
		}
				
		ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetRows);CHKERRQ(ierr);
		ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetColumns);CHKERRQ(ierr);
		
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, indexSetRows);CHKERRQ(ierr);
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, indexSetRows+1);CHKERRQ(ierr);
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, indexSetColumns);CHKERRQ(ierr);
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, indexSetColumns+1);CHKERRQ(ierr);
		
		for(int i=0; i<numSubMatrices; i++)
		{
			ierr = ISSort(indexSetRows[i]);CHKERRQ(ierr);
			ierr = ISSort(indexSetColumns[i]);CHKERRQ(ierr);
		}
		
		//}
		//if(rank == 0) numSubMatrices = p.size0; else numSubMatrices = 0;
		ierr = MatGetSubMatrices(inverseHamiltonian, numSubMatrices, indexSetRows, indexSetColumns, MAT_INITIAL_MATRIX, &localSubMatrices);CHKERRQ(ierr);

		
		//Fill gaml and gamrr
		
		for(int i=0; i<p.leadSize; i++)
		{
			ierr = MatGetRow(localSubMatrices[0], i, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int j=0; j<p.leadSize; j++)
				gaml[i][j] = rowValues[j];
			ierr = MatRestoreRow(localSubMatrices[0], i, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor		
		
		for(int j=0; j<p.leadSize; j++)
		{
			ierr = MatGetRow(localSubMatrices[1], j, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int i=0; i<p.leadSize; i++)
				gamrr[i][j] = conj(rowValues[i]);
			ierr = MatRestoreRow(localSubMatrices[1], j, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor
	
		//Destroy objects when finished
		//if(energy == energyMax)
		//{
		for(int i=0; i<numSubMatrices; i++)
		{
			ierr = ISDestroy(indexSetRows[i]);CHKERRQ(ierr);
			ierr = ISDestroy(indexSetColumns[i]);CHKERRQ(ierr);
			ierr = MatDestroy(localSubMatrices[i]);CHKERRQ(ierr);
		}
		ierr = PetscFree(indexSetRows);CHKERRQ(ierr);
		ierr = PetscFree(indexSetColumns);CHKERRQ(ierr);
		ierr = PetscFree(localSubMatrices);CHKERRQ(ierr);
		//} //endif
		
	} //endif
	else //serial
	{	
		for(int i=0; i<p.leadSize; i++)
		{
			ierr = MatGetRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int j=0; j<p.leadSize; j++)
				gaml[i][j] = rowValues[isr[j]-1];
			ierr = MatRestoreRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor
	
		for(int j=0; j<p.leadSize; j++)
		{
			ierr = MatGetRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int i=0; i<p.leadSize; i++)
				gamrr[i][j] = conj(rowValues[isr[i]-1]);
			ierr = MatRestoreRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor	
	} //endelse
	
		
	//Multiply the previous four leadSize x leadSize matrices together to get the transmissionFunction
	for(int i=0; i<p.leadSize; i++)
	{
		for (int j=0; j<p.leadSize; j++)
		{
			prodl[i][j] = complex<double>(0,0);
			prodr[i][j] = complex<double>(0,0);			
			
			for(int k=0; k<p.leadSize; k++)
			{
				prodl[i][j] += -gl[i][k].imag() * gaml[k][j];
				prodr[i][j] += -gr[i][k].imag() * gamrr[k][j];
			} //endfor
		} //endfor
	} //endfor	
	
	for(int i=0; i<p.leadSize; i++)
	{
		for(int j=0; j<p.leadSize; j++)
		{
			zres[i][j] = complex<double>(0,0);
			for(int k=0; k<p.leadSize; k++)
				zres[i][j] += prodl[i][k] * prodr[k][j];
		} //endfor
	} //endfor
	
	for(int i=0; i<p.leadSize; i++)
		transmissionFunction += zres[i][i];
	
	scaledTransmissionFunction = 4*pow(p.thop,4)*transmissionFunction;
	
	fermi1 = 1.0 / (exp((p.energy-emu1)/therm)+1.0);
	fermi2 = 1.0 / (exp((p.energy-emu2)/therm)+1.0);
	fermif = 2 * abs(transmissionFunction) * (fermi1-fermi2);
	
//	if(ranksg == 0)
//	{
		//Output to terminal
		ierr = PetscSynchronizedPrintf(PETSC_COMM_SUBGROUP, "%.5f \t%.5e\t   (%.5e, %.5e)\t%.5e\t%.0f\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real(), difftime(time(NULL), p.invertStart));CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_SUBGROUP);CHKERRQ(ierr);
	
		//Output to file
		ierr = PetscSynchronizedFPrintf(PETSC_COMM_SUBGROUP, outFile, "%.5f \t%.5e\t   (%.5e, %.5e)\t%.5e\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real());CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_SUBGROUP);CHKERRQ(ierr);
//	}
	
	//Send data to head node for later integration
	double IVvalue = 4 * pow(p.thop,4) * abs(transmissionFunction);
	
	if(size > 1)
	{
		if(rank != 0)
		{
			MPI_Send(&IVvalue, 1, MPI_DOUBLE, 0, currentStep, PETSC_COMM_WORLD);
		}
		else
		{
			MPI_Status MPIerr;
			IVvalues[currentStep] = IVvalue;
			int crossedZero = 0;
			for(int i=1; (i<size0n && currentStep+i<p.numEnergySteps); i++)
			{
				//Skip Recv if energy = 0
                if(abs((currentStep+i) * p.energyStep + p.energyMin) >= .9*p.energyStep)
				{	
					MPI_Recv(&IVvalues[currentStep+i-crossedZero], 1, MPI_DOUBLE, i, currentStep+i-crossedZero, PETSC_COMM_WORLD, &MPIerr);
				}
				else
					crossedZero = 1;
			} //endfor
		} //endelse
	} //endif
	else IVvalues[currentStep] = IVvalue;

	//ierr = MPI_Group_free(&PETSC_GROUP_WORLD);
	//ierr = MPI_Comm_free(&PETSC_COMM_GROUP);
	
	PetscFunctionReturn(0);
} //end RetrieveInformation


void permuteRow (complex<double> array[4][4], int rowA, int rowB)
{
	complex<double> tempVar;
	for(int j=0; j<4; j++)
	{
		tempVar = array[rowA][j];
		array[rowA][j] = array[rowB][j];
		array[rowB][j] = tempVar;
	}
} //end permuteRow


void permuteColumn (complex<double> array[4][4], int colA, int colB)
{
	complex<double> tempVar;
	for(int i=0; i<4; i++)
	{
		tempVar = array[i][colA];
		array[i][colA] = array[i][colB];
		array[i][colB] = tempVar;
	}
} //end permuteColumn

#endif /* INVERT_H_ */
