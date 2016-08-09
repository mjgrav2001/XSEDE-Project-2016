#ifndef INVERT_H_
#define INVERT_H_
#include <complex>
#include <iostream>
#include <ctime>
#include "petscksp.h"
using namespace std;

const double PI = 3.1415926535;

#undef __FUNCT__
#define __FUNCT__ "InvertHamiltonian"
PetscErrorCode InvertHamiltonian (Parameters& p, Mat hamiltonian, Mat inverseHamiltonian, Mat identityMat, PC pc, KSP ksp, 
									  block** blockHamiltonian)
{
	//PETSc Setup
	PetscErrorCode ierr;
		
	//Variable Declarations
	time(&p.invertStart);
	bool check;
	const PetscInt itor = p.ndiag/2;
    int ic[itor+1];
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
		
	//cout << " idim " << p.idim << " ndiag " << p.ndiag << endl;	
	
	//Fill Hamiltonian	
	for(int j=0; j<p.idim; j++)
	{
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				if(blockHamiltonian[j][1].getData(k, l) != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*j, blockHamiltonian[j][1].getData(k, l), 
									   INSERT_VALUES);CHKERRQ(ierr);
			} //endfor
		} //endfor
	} //endfor
	
	for(int j=0; j<p.idim-1; j++)
	{
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				if(blockHamiltonian[j][2].getData(k, l) != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k+p.ndiag*j, l+p.ndiag*(j-1), blockHamiltonian[j][2].getData(k, l), 
									   INSERT_VALUES);CHKERRQ(ierr);
				if(j == p.idim-2)
				{
					check = true;
					for(int i=1; i<=itor; i++)
					{
						if(check == true)
						{
							ic[i] = 2*(i-1);
							if(blockHamiltonian[j][2].getData(ic[i],ic[i]+1) != complex<double>(0,0))
								ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+2*(i-1), p.rows-2*p.ndiag+2*(i-1)+1, 
												   blockHamiltonian[j][2].getData(ic[i],ic[i]+1), INSERT_VALUES);CHKERRQ(ierr);
							check = !check;
						} //endif
						else
						{
							ic[i] = 2*i-1;
							if(blockHamiltonian[j][2].getData(ic[i],ic[i]-1) != complex<double>(0,0))
								ierr = MatSetValue(hamiltonian, p.rows-p.ndiag+2*(i-1)+1, p.rows-2*p.ndiag+2*(i-1), 
												   blockHamiltonian[j][2].getData(ic[i],ic[i]-1), INSERT_VALUES);CHKERRQ(ierr);
							check = !check;
						} //endelse
					} //endfor
				} //endif
			} //endfor
		} //endfor
	} //endfor
	
	for(int j=1; j<p.idim; j++)
	{
		for(int k=0; k<p.ndiag; k++)
		{
			for(int l=0; l<p.ndiag; l++)
			{
				if(blockHamiltonian[j][0].getData(k, l) != complex<double>(0,0))
					ierr = MatSetValue(hamiltonian, k+p.ndiag*(j-1), l+p.ndiag*j, blockHamiltonian[j][0].getData(k, l), 
									   INSERT_VALUES);CHKERRQ(ierr);
			} //endfor
		} //endfor
	} //endfor
	
	for(int k=0; k<p.ndiag; k++)
	{
		for (int l=0; l<p.ndiag; l++)
		{
			if(blockHamiltonian[0][0].getData(k, l) != complex<double>(0,0))
				ierr = MatSetValue(hamiltonian, k, l+p.ndiag*(p.idim-1), blockHamiltonian[0][0].getData(k, l), 
								   INSERT_VALUES);CHKERRQ(ierr);
			if(blockHamiltonian[p.idim-1][2].getData(k, l) != complex<double>(0,0))
				ierr = MatSetValue(hamiltonian, k+p.ndiag*(p.idim-1), l, blockHamiltonian[p.idim-1][2].getData(k, l), 
								   INSERT_VALUES);CHKERRQ(ierr);
		} //endfor
	} //endfor
	
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
	
	PetscFunctionReturn(0);
} //end InvertHamiltonian

#undef __FUNCT__
#define __FUNCT__ "RetrieveInformation"
PetscErrorCode RetrieveInformation (Parameters p, Mat inverseHamiltonian, double IVvalues[], FILE* outFile, int isl[], int isr[], complex<double> gl[4][4], complex<double> gr[4][4])
{
	//PETSc Setup
	const PetscScalar *rowValues;
	PetscScalar trace;
	PetscInt nonZeros, size, rank;
	PetscErrorCode ierr;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
	
	//Variable Declarations
	complex<double> gaml[4][4], gamrr[4][4];
	complex<double> prodl[4][4], prodr[4][4], zres[4][4];
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
	double densityOfStates = 0.0, therm = 0.03;
	double emu1 = -p.vsdMax/2.0, emu2 = p.vsdMax/2.0;
	int currentStep = (p.energy < 0.0 ? (p.energy - p.energyMin) / p.energyStep + 0.5 : (p.energy - p.energyMin) / p.energyStep - 0.5);
	
	ierr = MatGetTrace(inverseHamiltonian, &trace);CHKERRQ(ierr);
	densityOfStates = -1.0 * trace.imag() / PI;
	
	for(int i=0; i<p.leadSize; i++)
	{
		ierr = MatGetRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		
//		cout << " i=" << i << " isl[i]=" << isl[i] << endl;
		
		for(int j=0; j<p.leadSize; j++)
		{
			gaml[i][j] = rowValues[isr[j]-1];
//			cout << " j=" << j << " isr[j]=" << isr[j] << endl;
			
		}
		ierr = MatRestoreRow(inverseHamiltonian, isl[i]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
	} //endfor
	
	for(int j=0; j<p.leadSize; j++)
	{
		ierr = MatGetRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
		for(int i=0; i<p.leadSize; i++)
			gamrr[i][j] = conj(rowValues[isr[i]-1]);
		ierr = MatRestoreRow(inverseHamiltonian, isl[j]-1, &nonZeros, PETSC_NULL, &rowValues);CHKERRQ(ierr);
	} //endfor

//	
//	for(int i=0; i<p.leadSize; i++)
//	{
//		for(int j=0; j<p.leadSize; j++)
//		{
//			cout << ' ' << i+1 <<  ' ' << j+1 <<  ' ' << gaml[i][j] << ' ' << gamrr[i][j] << endl;
//		}
//	} //endfor
//
	
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
	
	scaledTransmissionFunction = 4.0*pow(p.thop,4)*transmissionFunction;
//	scaledTransmissionFunction = 4.0*pow(-0.25,4)*transmissionFunction;
	
	fermi1 = 1.0 / (exp((p.energy-emu1)/therm)+1.0);
	fermi2 = 1.0 / (exp((p.energy-emu2)/therm)+1.0);
	fermif = 2.0 * abs(transmissionFunction) * (fermi1-fermi2);
	
	//Output to terminal
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%.5f \t%.5e\t (%.5e, %.5e)\t%.5e\t%.0f\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real(), difftime(time(NULL), p.invertStart));CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	//Output to file
	ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "%.5f \t%.5e\t (%.5e, %.5e)\t%.5e\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real());CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	
	//Send data to head node for later integration
	double IVvalue = 4.0 * pow(p.thop,4) * abs(transmissionFunction);
//	double IVvalue = 4.0 * pow(-0.25,4) * abs(transmissionFunction);
	if(size > 1)
	{
		if(rank != 0)
			MPI_Send(&IVvalue, 1, MPI_DOUBLE, 0, currentStep, PETSC_COMM_WORLD);
		else
		{
			MPI_Status MPIerr;
			IVvalues[currentStep] = IVvalue;
			int crossedZero = 0;
			for(int i=1; (i<size && currentStep+i<p.numEnergySteps); i++)
			{
				//Skip Recv if energy = 0
                if(abs((currentStep+i) * p.energyStep + p.energyMin) >= 0.9*p.energyStep)
					MPI_Recv(&IVvalues[currentStep+i-crossedZero], 1, MPI_DOUBLE, i, currentStep+i-crossedZero, PETSC_COMM_WORLD, &MPIerr);
				else
					crossedZero = 1;
			} //endfor
		} //endelse
	} //endif
	else IVvalues[currentStep] = IVvalue;
	
	PetscFunctionReturn(0);
} //end RetrieveInformation

#endif /* INVERT_H_ */