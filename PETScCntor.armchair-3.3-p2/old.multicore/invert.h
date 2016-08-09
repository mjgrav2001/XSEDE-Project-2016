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

#undef __FUNCT__
#define __FUNCT__ "InvertHamiltonian"
PetscErrorCode InvertHamiltonian (Parameters& p, Mat hamiltonian, Mat inverseHamiltonian, Mat identityMat, PC pc, KSP ksp, 
								  int isl180[], int isr180[], int isl[], int isr[], complex<double> gl[4][4], complex<double> gr[4][4], 
								  double atomCoordinates[][3], double disCoordinates[][3], block blockHamiltonian)
{
	//PETSc Setup
	PetscErrorCode ierr;
	
	//Variable Declarations';
	time(&p.invertStart);
	bool check;
	const PetscInt itor = p.idim/2;
    int ici,ici1;
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
	complex<double> ham;
	
	//Fill Hamiltonian	
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
	
	
	for(int j=1; j<p.ndiag; j++)
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
	
	} //endfor

	
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
	PetscInt nonZeros, sizen, rankn, size, rank;
	PetscErrorCode ierr;

	//Variable Declarations
	complex<double> gaml[4][4], gamrr[4][4];
	complex<double> prodl[4][4], prodr[4][4], zres[4][4];
	complex<double> fermi1, fermi2, fermif, transmissionFunction;
	complex<double> scaledTransmissionFunction;
	double densityOfStates = 0.0, therm = .03;
	double emu1 = -p.vsdMax/2.0, emu2 = p.vsdMax/2.0;
	int currentStep = (p.energy < 0.0 ? (p.energy - p.energyMin) / p.energyStep + .5 : (p.energy - p.energyMin) / p.energyStep - .5);

	//******************************************************************************************************
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
	//******************************************************************************************************
	
	int n;
	n =  (int) size / p.size0;
	
	int k;
	int ranksn[n];
	for (k=0; k<n; k++)
	{
		ranksn[k] = k*p.size0;
	}

	//cout <<  " rank " << rank << " size " << size << endl;
	//******************************************************************************************************
	MPI_Group PETSC_GROUP_WORLD, PETSC_GROUP;
	MPI_Comm PETSC_COMM_GROUP;
	
	MPI_Comm_group(PETSC_COMM_WORLD, &PETSC_GROUP_WORLD);
	MPI_Group_incl(PETSC_GROUP_WORLD, n, ranksn, &PETSC_GROUP); /* local */
	
	MPI_Comm_create(PETSC_COMM_WORLD, PETSC_GROUP_WORLD, &PETSC_COMM_GROUP);
	
	ierr = MPI_Comm_size(PETSC_COMM_GROUP, &sizen);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_GROUP, &rankn);CHKERRQ(ierr);
	//******************************************************************************************************
	//cout <<  " rankn " << rankn << " sizen " << sizen << endl;
	
	ierr = MatGetTrace(inverseHamiltonian, &trace);CHKERRQ(ierr);
	densityOfStates = -1.0 * trace.imag() / PI;
	
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
	
	//Output to terminal
    ierr = PetscSynchronizedPrintf(PETSC_COMM_GROUP, "%.5f \t%.5e\t   (%.5e, %.5e)\t%.5e\t%.0f\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real(), difftime(time(NULL), p.invertStart));CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_GROUP);CHKERRQ(ierr);
	
	//Output to file
	ierr = PetscSynchronizedFPrintf(PETSC_COMM_GROUP, outFile, "%.5f \t%.5e\t   (%.5e, %.5e)\t%.5e\n", p.energy, densityOfStates, scaledTransmissionFunction.real(), scaledTransmissionFunction.imag(), fermif.real());CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_GROUP);CHKERRQ(ierr);
	
	//Send data to head node for later integration
	double IVvalue = 4 * pow(p.thop,4) * abs(transmissionFunction);
	if(sizen > 1)
	{
		if(rankn != 0)
			MPI_Send(&IVvalue, 1, MPI_DOUBLE, 0, currentStep, PETSC_COMM_GROUP);
		else
		{
			MPI_Status MPIerr;
			IVvalues[currentStep] = IVvalue;
			int crossedZero = 0;
			for(int i=1; (i<sizen && currentStep+i<p.numEnergySteps); i++)
			{
				//Skip Recv if energy = 0
                if(abs((currentStep+i) * p.energyStep + p.energyMin) >= .9*p.energyStep)
					MPI_Recv(&IVvalues[currentStep+i-crossedZero], 1, MPI_DOUBLE, i, currentStep+i-crossedZero, PETSC_COMM_GROUP, &MPIerr);
				else
					crossedZero = 1;
			} //endfor
		} //endelse
	} //endif
	else IVvalues[currentStep] = IVvalue;

	PetscFunctionReturn(0);
} //end RetrieveInformation

#endif /* INVERT_H_ */
