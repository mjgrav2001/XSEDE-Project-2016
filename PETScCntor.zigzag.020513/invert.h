#ifndef INVERT_H_
#define INVERT_H_
#include "petscksp.h"
#include "petsc.h"
#include "classes.h"
#include "general.h"
#include "grid_jobs.h"
#include <iostream>
#include <complex>
#include <new>
#include <ctime>
#include <fstream>

using namespace std;

const double PI = 3.1415926535;

#undef __FUNCT__
#define __FUNCT__ "InvertHamiltonian"
PetscErrorCode InvertHamiltonian (MPI_Comm PETSC_COMM_WORLD, Parameters& p, Mat hamiltonian, Mat inverseHamiltonian, Mat identityMat, PC pc, KSP ksp, int *isl180, int *isr180, int *isl, int *isr, complex<double> **gl, complex<double> **gr, double atomCoordinates[][3], double disCoordinates[][3])
{
	//Declarations
	time_t programStart;
	time(&programStart);
	
	//PETSc Setup
	PetscErrorCode ierr;
	
	//Variable Declarations';
	time(&p.invertStart);
	bool check;
	PetscInt ndiag1 = p.ndiag;
	PetscInt itor = p.ndiag/2;
	
	block blockHamiltonian;
    
	int ici,ici1;
	complex<double> fermi1, fermi2, fermif;
	PetscScalar *ham, *ham0;	
	ierr = PetscMalloc(itor*itor*sizeof(PetscScalar), &ham0);CHKERRQ(ierr);		
	ierr = PetscMalloc(ndiag1*ndiag1*sizeof(PetscScalar), &ham);CHKERRQ(ierr);
	
	int m = (int) p.idim / p.myGroupSize;
	
	for(int j=p.myGroupRank*m; j<(p.myGroupRank+1)*m; j++)
	{
			blockHamiltonian.set_zeros(p);
			blockHamiltonian.set_values(p, j, 1, 'D', atomCoordinates, disCoordinates);
			blockHamiltonian.AttachLeads(p, j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);			
			for(int k=0; k<p.ndiag; k++)
			{
				for(int l=0; l<p.ndiag; l++)
				{
					ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
					if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
						ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*j, ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
					}
				} //endfor
			} //endfor
			//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
			//							j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);				
	}
	
	if(p.idim % p.myGroupSize != 0)
	{
		if(p.myGroupRank == 0)
		{
			for(int j=m*p.myGroupSize; j<p.idim; j++)
			{
					blockHamiltonian.set_zeros(p);
					blockHamiltonian.set_values(p, j, 1, 'D', atomCoordinates, disCoordinates);
					blockHamiltonian.AttachLeads(p, j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);			
					for(int k=0; k<p.ndiag; k++)
					{
						for(int l=0; l<p.ndiag; l++)
						{
							ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);	
							if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
								ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*j, ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
							}
						} //endfor
					} //endfor
					//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
					//							j, 1, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
			}
		}
	}
	
	//******************************************************************************************************
	for(int j=p.myGroupRank*m; j<(p.myGroupRank+1)*m; j++)
	{
		if(j<p.idim-1)
		{
			blockHamiltonian.set_zeros(p);
			blockHamiltonian.set_values(p, j, 2, 'W', atomCoordinates, disCoordinates);
			blockHamiltonian.AttachLeads(p, j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
			for(int k=0; k<p.ndiag; k++)
			{
				for(int l=0; l<p.ndiag; l++)
				{
					ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
					if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
						ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*(j-1), ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
					}
					
					if(j == p.idim-2)
					{
						check = true;
						for(int i=1; i<=itor; i++)
						{
							if(check == true)
							{
								ici = 2*(i-1);
								ici1 = ici+1;
								ham0[(i-1)*itor+i-1] = blockHamiltonian.get_data(ici,ici1);
								if(ham0[(i-1)*itor+i-1]!= complex<double>(0,0)) {
									ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+ici, p.rows-2*p.ndiag+ici1, ham0[(i-1)*itor+i-1], INSERT_VALUES);CHKERRQ(ierr);
								}
								check = !check;
							} //endif
							else
							{
								ici = 2*i-1;
								ici1 = ici-1;
								ham0[(i-1)*itor+i-1] = blockHamiltonian.get_data(ici,ici1);
								if(ham0[(i-1)*itor+i-1] != complex<double>(0,0)) {
									ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+ici, p.rows-2*p.ndiag+ici1, ham0[(i-1)*itor+i-1], INSERT_VALUES);CHKERRQ(ierr);
								}
								check = !check;
							} //endelse
						} //endfor
					} //endif
				} //endfor
			} //endfor
			//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
			//							 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);				
		} //endfor
	}

	if(p.idim % p.myGroupSize != 0)
	{
		if(p.myGroupRank == 0)
		{
			for(int j=m*p.myGroupSize; j<p.idim-1; j++)
			{
				blockHamiltonian.set_zeros(p);
				blockHamiltonian.set_values(p, j, 2, 'W', atomCoordinates, disCoordinates);
				blockHamiltonian.AttachLeads(p, j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
				for(int k=0; k<p.ndiag; k++)
				{
					for(int l=0; l<p.ndiag; l++)
					{
						ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
						if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
							ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*(j-1), ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
						}
						if(j == p.idim-2)
						{
							check = true;
							for(int i=1; i<=itor; i++)
							{
								if(check == true)
								{
									ici = 2*(i-1);
									ici1 = ici+1;
									ham0[(i-1)*itor+i-1] = blockHamiltonian.get_data(ici,ici1);
									if(ham0[(i-1)*itor+i-1]!= complex<double>(0,0)) {
										ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+ici, p.rows-2*p.ndiag+ici1, ham0[(i-1)*itor+i-1], INSERT_VALUES);CHKERRQ(ierr);
									}
									check = !check;
								} //endif
								else
								{
									ici = 2*i-1;
									ici1 = ici-1;
									ham0[(i-1)*itor+i-1] = blockHamiltonian.get_data(ici,ici1);
									if(ham0[(i-1)*itor+i-1] != complex<double>(0,0)) {
										ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+ici, p.rows-2*p.ndiag+ici1, ham0[(i-1)*itor+i-1], INSERT_VALUES);CHKERRQ(ierr);
									}
									check = !check;										
								} //endelse
							} //endfor
						} //endif
					} //endfor
				} //endfor
				//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
				//							 j, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
			}
		}
	}
	
	//******************************************************************************************************
	for(int j=p.myGroupRank*m; j<(p.myGroupRank+1)*m; j++)
	{
		blockHamiltonian.set_zeros(p);
		blockHamiltonian.set_values(p, j-1, 0, 'T', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p, j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);	
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
				if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
					ierr = MatSetValue(hamiltonian, k+p.ndiag*(j-1), l+p.ndiag*j, ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
				}
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);			
	}
	
	if(p.idim % p.myGroupSize != 0)
	{
		if(p.myGroupRank == 0)
		{
			for(int j=m*p.myGroupSize; j<p.idim; j++)
			{
				blockHamiltonian.set_zeros(p);
				blockHamiltonian.set_values(p, j-1, 0, 'T', atomCoordinates, disCoordinates);
				blockHamiltonian.AttachLeads(p, j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);	
				for(int k=0; k<p.ndiag; k++)
				{
					for(int l=0; l<p.ndiag; l++)
					{
						ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
						if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
							ierr = MatSetValue(hamiltonian, k+p.ndiag*(j-1), l+p.ndiag*j, ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
						}
					} //endfor
				} //endfor
				//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
				//							 j-1, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 1);		
			}
		}
	}

	//******************************************************************************************************
	if(p.myGroupRank == 0)
	{
		blockHamiltonian.set_zeros(p);
		blockHamiltonian.set_values(p, 0, 0, 'U', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p, 0, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);		
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
				if(ham[k*p.ndiag+l]!= complex<double>(0,0)) {
					ierr = MatSetValue(hamiltonian, k, l+p.ndiag*(p.idim-1), ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
				}
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 0, 0, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	

		blockHamiltonian.set_zeros(p);
		blockHamiltonian.set_values(p, p.idim-1, 2, 'Q', atomCoordinates, disCoordinates);
		blockHamiltonian.AttachLeads(p, p.idim-1, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);		
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				ham[k*p.ndiag+l] = blockHamiltonian.get_data(k,l);
				if(ham[k*p.ndiag+l] != complex<double>(0,0)) {
					ierr = MatSetValue(hamiltonian, k+p.ndiag*(p.idim-1), l, ham[k*p.ndiag+l], INSERT_VALUES);CHKERRQ(ierr);
				}
			} //endfor
		} //endfor
		//blockHamiltonian.DetachLeads(p.leadSize, p.rows, p.ndiag, p.idim, p.vhop, p.energy, p.BField, p.thop, blockHamiltonian, 
		//							 p.idim-1, 2, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates, 0);	
	}
	
	//******************************************************************************************************
	//Assemble PETSc Hamiltonian
	ierr = MatAssemblyBegin(hamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(hamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Reset PETSc Solver
	ierr = KSPSetOperators(ksp, hamiltonian, hamiltonian, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);
	
	//Factor hamiltonian
	ierr = PCFactorGetMatrix(pc, &hamiltonian);CHKERRQ(ierr);
	
	//Solve hamiltonian * inverseHamiltonian = identityMat for inverseHamiltonian
	ierr = MatMatSolve(hamiltonian, identityMat, inverseHamiltonian);CHKERRQ(ierr);
	
	//******************************************************************************************************
	//Print run time
	//if(p.myWorldRank == 0)
	//{
	//	int programTime = (int) difftime(time(NULL), programStart);
	//	int seconds = programTime % 60;
	//	int minutes = (programTime / 60) % 60;
	//	int hours = programTime / 3600;
	//	cout << "\nTotal run time: ";
	//	if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
	//	else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
	//	else cout << seconds << " seconds" << endl;
	//}
	
	ierr = PetscFree(ham);CHKERRQ(ierr);
	ierr = PetscFree(ham0);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
} //end InvertHamiltonian


//******************************************************************************************************
#undef __FUNCT__
#define __FUNCT__ "RetrieveInformation"
PetscErrorCode RetrieveInformation (Parameters& p, int& numSubMatrices, IS *indexSetRows, IS *indexSetColumns,
									Mat inverseHamiltonian, int *isl, int *isr, double& densityOfStates, complex<double> **gl, complex<double> **gr, 
									complex<double> **gaml, complex<double> **gamrr)
{
	//Declarations
	time_t programStart;
	time(&programStart);
	
	//PETSc Setup
	static Mat *localSubMatrices;
	const PetscScalar *rowValues;
	
	PetscScalar trace;
	PetscInt nonZeros;
	PetscErrorCode ierr;

	ierr = MatGetTrace(inverseHamiltonian, &trace);CHKERRQ(ierr);
	densityOfStates = -1.0 * trace.imag() / PI;
	
	//******************************************************************************************************
	//if(p.myGroupSize > 1 || p.worldSize > 1) //parallel
	if(p.myGroupSize > 1) //parallel
	{	
		int lsize=p.leadSize;	
		ierr = PetscMalloc(lsize*sizeof(PetscScalar), &rowValues);CHKERRQ(ierr);
		
		ierr = MatGetSubMatrices(inverseHamiltonian, numSubMatrices, indexSetRows, indexSetColumns, MAT_INITIAL_MATRIX, &localSubMatrices);CHKERRQ(ierr);
		
		//Fill gaml and gamrr
		for(int i=0; i<p.leadSize; i++) {
			ierr = MatGetRow(localSubMatrices[0], i, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int j=0; j<p.leadSize; j++) {
				gaml[i][j] = rowValues[j];
				//cout <<  "gaml[" << i << "][" << j << "] = " << gaml[i][j] << endl;
			}
			ierr = MatRestoreRow(localSubMatrices[0], i, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		}
		
		for(int j=0; j<p.leadSize; j++) {
			ierr = MatGetRow(localSubMatrices[0], j, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int i=0; i<p.leadSize; i++) {
				gamrr[i][j] = conj(rowValues[i]);
				//cout <<  "gamrr[" << i << "][" << j << "] = " << gamrr[i][j] << endl;
			}
			ierr = MatRestoreRow(localSubMatrices[0], j, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		}
		
		//******************************************************************************************************	
		//Determine permutation needed
			
		//allocate memory for elements of each column
		int *islPerm, *isrPerm;
		islPerm = new int[lsize];
		isrPerm = new int[lsize];
		
		for(int i=0; i<lsize; i++) {
			isrPerm[i] = 0;
			islPerm[i] = 0;
			
			for(int j=0; j<lsize; j++) {
				if(isr[j] < isr[i])	isrPerm[i]++;
				if(isl[j] < isl[i])	islPerm[i]++;
			}
		}
		
		//Permute rows (islPerm) and columns (isrPerm)
		for(int i=0; i<lsize; i++) {
			int tempVar;
			if(islPerm[i] != i) {
				permuteRow(lsize, gaml, i, islPerm[i]);
				permuteColumn(lsize, gamrr, i, islPerm[i]);
				tempVar = islPerm[i];
				islPerm[i] = i;
				islPerm[tempVar] = tempVar;
			}
		}
		
		for(int j=0; j<lsize; j++) {
			int tempVar;			
			if(isrPerm[j] != j) {
				permuteColumn(lsize, gaml, j, isrPerm[j]);				
				permuteRow(lsize, gamrr, j, isrPerm[j]);
				tempVar = isrPerm[j];
				isrPerm[j] = j;
				isrPerm[tempVar] = tempVar;
			}
			
			//if(p.myWorldRank == 0) {
				//for(int i=0; i<p.leadSize; i++) {		
					//cout <<  "gaml[" << i << "][" << j << "] = " << gaml[i][j] << endl;
					//cout <<  "gamrr[" << i << "][" << j << "] = " << gamrr[i][j] << endl;
				//}
			//}
		}
		
		//de-allocate memory
		delete[] islPerm;
		delete[] isrPerm;
		
		ierr = PetscFree(rowValues);CHKERRQ(ierr);
	}
	else { //serial
		int lsize = p.columns;
		ierr = PetscMalloc(lsize*sizeof(PetscScalar), &rowValues);CHKERRQ(ierr);
		
		for(int i=0; i<p.leadSize; i++) {
			ierr = MatGetRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int j=0; j<p.leadSize; j++) {
				gaml[i][j] = rowValues[isr[j]-1];
				//cout <<  "gaml[" << i << "][" << j << "] = " << gaml[i][j] << endl;
			}
			ierr = MatRestoreRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor
		
		for(int j=0; j<p.leadSize; j++) {
			ierr = MatGetRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
			for(int i=0; i<p.leadSize; i++) {
				gamrr[i][j] = conj(rowValues[isr[i]-1]);
				//cout <<  "gamrr[" << i << "][" << j << "] = " << gamrr[i][j] << endl;
			}
			ierr = MatRestoreRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		} //endfor
		
		ierr = PetscFree(rowValues);CHKERRQ(ierr);	
	}

	//cout <<  "3. p.energy " << p.energy << " p.myWorldRank " << p.myWorldRank << endl;

	//Destroy PETSc Objects and free memory
	ierr = MatDestroy(localSubMatrices[0]);CHKERRQ(ierr);
 	ierr = PetscFree(localSubMatrices);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
} //end RetrieveInformation


//******************************************************************************************************
#undef __FUNCT__
#define __FUNCT__ "CalculateTransport"
PetscErrorCode CalculateTransport (MPI_Comm PETSC_COMM_WORLD, Parameters& p, double& densityOfStates, complex<double> **gl, complex<double> **gr, complex<double> **gaml, complex<double> **gamrr, double IVvalues[], FILE* outFile)
{
	//Declarations
	time_t programStart;
	time(&programStart);
	
	//PETSc setup	
	PetscErrorCode ierr;
		
	//Variable Declarations
	double therm = .03;
	double emu1 = -p.vsdMax/2.0, emu2 = p.vsdMax/2.0;
	complex<double> transmissionFunction;
	double IVvalue;
	
	//******************************************************************************************************
	//(step 0) declare arrays
	complex<double> **prodl, **prodr, **zres;
	
	//******************************************************************************************************
	//(step 1) allocate memory for array of elements of each column of coupling matrices
	prodl = new complex<double>*[p.leadSize];
	prodr = new complex<double>*[p.leadSize];
	zres = new complex<double>*[p.leadSize];
	
	//(step 2) allocate memory for array of elements of each row of coupling matrices
	for( int i = 0; i < p.leadSize; i++) {
		*(prodl + i) = new complex<double>[p.leadSize];
		*(prodr + i) = new complex<double>[p.leadSize];
		*(zres + i) = new complex<double>[p.leadSize];
	}
	
	int currentStep = (p.energy < 0.0 ? (p.energy - p.energyMin) / p.energyStep + .5 : (p.energy - p.energyMin) / p.energyStep - .5);
		
	//******************************************************************************************************	
	//Multiply the previous four leadSize x leadSize matrices together to get the transmissionFunction
	for(int i=0; i<p.leadSize; i++)
	{
		for(int j=0; j<p.leadSize; j++)
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

	complex<double> scaledTransmissionFunction = 4.0*pow(p.thop,4)*transmissionFunction;	
	complex<double> fermi1 = 1.0 / (exp((p.energy-emu1)/therm)+1.0);
	complex<double> fermi2 = 1.0 / (exp((p.energy-emu2)/therm)+1.0);
	complex<double> fermif = 2.0 * abs(transmissionFunction) * (fermi1-fermi2);
	
	IVvalue = 4.0 * pow(p.thop,4) * abs(transmissionFunction);
	
	//******************************************************************************************************	
	//De-allocatte memory
	for(int i = 0; i < p.leadSize; i++)
	{
		delete[] *(prodl + i);
		delete[] *(prodr + i);
		delete[] *(zres + i);
	}
	delete[] prodl;
	delete[] prodr;
	delete[] zres;

	//int programTime;
	//int seconds;
	//int minutes;
	//int hours;
	
	//******************************************************************************************************
	//Create output on screen and in output file
	if((p.myWorldRank % p.groupSize) == 0)
	{	
		//Output to terminal
		ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%.5f \t%.5e\t (%.5e, %.5e)\t%.5e\t%.0f\n", p.energy, densityOfStates, 
									   scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real(), 
									   difftime(time(NULL), p.invertStart));CHKERRQ(ierr);
		
		//programTime = (int) difftime(time(NULL), programStart);
		//seconds = programTime % 60;
		//minutes = (programTime / 60) % 60;
		//hours = programTime / 3600;
		//cout << "\nTotal run time: ";
		//if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
		//else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
		//else cout << seconds << " seconds" << endl;
		
	}		
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	
	if((p.myWorldRank % p.groupSize) == 0)
	{
		//Output to file
		ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "%.5f \t%.5e\t (%.5e, %.5e)\t%.5e\n", p.energy, densityOfStates, 
										scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real());CHKERRQ(ierr);
		
		//programTime = (int) difftime(time(NULL), programStart);
		//seconds = programTime % 60;
		//minutes = (programTime / 60) % 60;
		//hours = programTime / 3600;
		//cout << "\nTotal run time: ";
		//if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
		//else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
		//else cout << seconds << " seconds" << endl;		
		
	}		
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	//cout <<  "3. p.energy " << p.energy << " p.myWorldRank " << p.myWorldRank << endl;
	
	//******************************************************************************************************
	//Send data to head node for later integration
	MPI_Status MPIerr;
	int crossedZero = 0;
	
	if(p.procsInUse > 1)
	{
		if(p.myWorldRank != 0)
		{
			if((p.myWorldRank % p.groupSize) == 0)
			{
				MPI_Send(&IVvalue, 1, MPI_DOUBLE, 0, currentStep, PETSC_COMM_WORLD);
			}
		}
		else
		{
			IVvalues[currentStep] = IVvalue;
			for(int i=1; (i < p.numGroups && currentStep+i*p.groupSize < p.numEnergySteps); i++)
			{
				//Skip Recv if energy = 0
			  if(abs((currentStep+i*p.groupSize) * p.energyStep + p.energyMin) >= .9*p.energyStep)
				{	
					
				  //programTime = (int) difftime(time(NULL), programStart);
				  //seconds = programTime % 60;
				  //minutes = (programTime / 60) % 60;
				  //hours = programTime / 3600;
				  //cout << "\nTotal run time: ";
				  //if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
				  //else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
				  //else cout << seconds << " seconds" << endl;
					
					MPI_Recv(&IVvalues[currentStep+i-crossedZero], 1, MPI_DOUBLE, i*p.groupSize, 
							 currentStep+i-crossedZero, PETSC_COMM_WORLD, &MPIerr);
					
					//programTime = (int) difftime(time(NULL), programStart);
					//seconds = programTime % 60;
					//minutes = (programTime / 60) % 60;
					//hours = programTime / 3600;
					//cout << "\nTotal run time: ";
					//if(hours != 0) cout << hours << " hours " << minutes << " minutes " << seconds << " seconds" << endl;
					//else if(minutes != 0) cout << minutes << " minutes " << seconds << " seconds" << endl;
					//else cout << seconds << " seconds" << endl;					
					
				}
				else
					crossedZero = 1;
			} //endfor
		} //endelse
	} //endif
	else IVvalues[currentStep] = IVvalue;
	
	PetscFunctionReturn(0);
} //end CalculateTransport

#endif /* INVERT_H_ */
