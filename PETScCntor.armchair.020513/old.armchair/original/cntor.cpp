#include "petscksp.h"
#include "classes.h"
#include "grid_jobs.h"
#include "invert.h"
#include "general.h"
#include "integration.h"
#include <iostream>
#include <complex>
#include <new>
#include <ctime>
#include <fstream>
using namespace std;

int main(int argc,char **args)
{
	//Declarations
	time_t programStart;
	time(&programStart);
	PetscInitialize(&argc, &args);
	PC pc;
	KSP ksp;
	PetscInt size, rank;
	PetscErrorCode ierr;
	Mat hamiltonian, inverseHamiltonian, identityMat;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
	complex<double> gl[4][4], gr[4][4], dos2;
	double corrsd, vsd;
	
	//Numerical accuracies
	double eta1=1e-18;
	
	//Crtical values:
	double bnatural = 0.004762;
	//	Bfield = bnatural / 16.;
	
	//Setting all parameters in a struct for easier passing to functions and a central location to make changes
	//IMPORTANT: Subroutines will have their own *copy* of this class unless passed by reference.
	//			 If passed by reference, changes will be seen by all subsequent uses of the class.
	Parameters p(      /*rows*/ 4000,
				    /*columns*/ 4000,
				       /*idim*/ 40,
			      /*energyMin*/ -0.2,
				  /*energyMax*/ 0.2,
				 /*energyStep*/ 0.00005,
				       /*thop*/ -0.25,
				     /*vsdMin*/ 0.0,
				     /*vsdMax*/ 0.1,
				   /*vsdSteps*/ 100,
				       /*vhop*/ -3.1,
				     /*BField*/ 0.0,
				   /*leadSize*/ 4,
				       /*mtop*/ 100,
				       /*ntop*/ 100,
				          /*a*/ 2.0);
	
	const int leadSize = 4;
	double IVvalues[p.numEnergySteps];
	
	cout <<  " Bfield " << p.BField << endl;
	cout << " rows " << p.rows << " idim " << p.idim << " ndiag " << p.ndiag << " numEnergySteps " << p.numEnergySteps << endl;
	
	//block blockHamiltonian[p.ndiag][3];
	block** blockHamiltonian;
	//(step 1) allocate memory for array of elements of column
	blockHamiltonian = (block** ) calloc(p.ndiag, sizeof(block* ));
	//(step 2) allocate memory for array of elements of each row
	for( int i = 0; i < p.ndiag; i++)
		*(blockHamiltonian + i) = (block* ) calloc(3, sizeof(block));
	
	//double atomCoordinates[p.rows][3];	
	double **atomCoordinates;	
	//(step 1) allocate memory for array of elements of column
	atomCoordinates = (double** ) calloc(p.rows, sizeof(double* ));
	//(step 2) allocate memory for array of elements of each row
	for( int i = 0; i < p.rows; i++)
		*(atomCoordinates + i) = (double* ) calloc(3, sizeof(double));
	
	//double disCoordinates[p.rows][3];	
	double **disCoordinates;	
	//(step 1) allocate memory for array of elements of column
	disCoordinates = (double** )calloc(p.rows, sizeof(double* ));
	//(step 2) allocate memory for array of elements of each row
	for( int i = 0; i < p.rows; i++)
		*(disCoordinates + i) = (double* )calloc(3, sizeof(double));
	
	//double p.BField = (8./32.)*(0.0255664/2.);

	//double p.BField = (8./32.)*(0.0255664/2.)
	int	totalSteps = ((p.energyMin <= 0 && p.energyMax >= 0) ? p.numEnergySteps + 1 : p.numEnergySteps),
		isl[leadSize] = {1, 2, p.imax1, p.imax2},
		isr[leadSize] = {p.imaxd21, p.imaxd22, p.imaxd21i, p.imaxd22i},
		isl180[leadSize] = {1, 2, p.imax1, p.imax2},
		isr180[leadSize] = {p.imaxd21, p.imaxd22, p.imaxd21i, p.imaxd22i};
		
	//Open output files
	FILE *outFile, *IVoutput;
	bool failFlag = false;
	
	ierr = PetscFOpen(PETSC_COMM_WORLD, "DETEfEoutput.txt", "w", &outFile);CHKERRQ(ierr);
	ierr = PetscFOpen(PETSC_COMM_WORLD,"IVoutput.txt", "w", &IVoutput);CHKERRQ(ierr);
	if(rank == 0 && (outFile == NULL || IVoutput == NULL))
	{
		cout << "Error: Output files not opened!\n";
		failFlag = true;
	} //endif
	if(failFlag)
	{
		ierr = PetscFinalize();CHKERRQ(ierr);
		return 0;
	}
	
	//Print start time
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\nStart time: %s", ctime(&programStart));CHKERRQ(ierr);
	
	//Checking number of processors
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Processors in use: %d\n", size);CHKERRQ(ierr);
	
//******************************************************************************************************
	//Setup hamiltonian
	ierr = MatCreate(PETSC_COMM_SELF, &hamiltonian);CHKERRQ(ierr);
	ierr = MatSetType(hamiltonian, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(hamiltonian, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr);
	
	//Setup PETSc Identity
	ierr = MatCreate(PETSC_COMM_SELF, &identityMat);CHKERRQ(ierr);
	ierr = MatSetSizes(identityMat, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr);
	ierr = MatSetType(identityMat, MATDENSE);CHKERRQ(ierr);
	for(int i=0; i<p.rows; i++)
		ierr = MatSetValue(identityMat, i, i, 1.0, INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//Setup PETSc inverseMat
	ierr = MatDuplicate(identityMat, MAT_DO_NOT_COPY_VALUES, &inverseHamiltonian);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//Setup PETSc Solver
	ierr = KSPCreate(PETSC_COMM_SELF, &ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, hamiltonian, hamiltonian, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
	ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	ierr = PCFactorSetZeroPivot(pc, eta1);CHKERRQ(ierr);
//******************************************************************************************************
	
	//Read in and setup system
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Energy range: %.5f eV to %.5f eV\n", p.energyMin, p.energyMax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "BField = %f\n", p.BField);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\nEnergy\t\tDensityofStates\t   ScaledTransmissionFunction\tRe(fermif)\tTime\n\n");CHKERRQ(ierr);	
	
	ReadInCoordinates(p.rows, atomCoordinates);	
	ReadInDisplaced(p.rows, disCoordinates);
	
	InitGridOffDiag(p, blockHamiltonian, atomCoordinates, disCoordinates);		

	//Loop over energy values
	for(p.energy=p.energyMin+rank*p.energyStep; p.energy<(p.energyMax+.9*p.energyStep); p.energy+=p.energyStep*size)
	{
        //Skipping energy = 0
		if(abs(p.energy) >= .9*p.energyStep)
		{
			InitGridDiag(p, blockHamiltonian, atomCoordinates, disCoordinates);
			wiregf(p, isl180, isr180, isl, isr, gl, gr, disCoordinates);
			AttachLeads(p, blockHamiltonian, isl180, isr180, isl, isr, gl, gr);
			ierr = InvertHamiltonian(p, hamiltonian, inverseHamiltonian, identityMat, pc, ksp, blockHamiltonian);CHKERRQ(ierr);
			ierr = RetrieveInformation(p, inverseHamiltonian, IVvalues, outFile, isl, isr, gl, gr);CHKERRQ(ierr);
			DetachLeads(p, blockHamiltonian, isl180, isr180, isl, isr, gl, gr);
		} //endif
		else
		{
            //Empty sync when skipping energy = 0 to ensure output order
            ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
            ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
            ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
            ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	    }
	} //endfor
	ierr = PetscFClose(PETSC_COMM_WORLD, outFile);CHKERRQ(ierr);
	
	//Determine if extra step is needed for synchronized prints
	if(totalSteps % size != 0 && rank >= totalSteps % size)
	{
        //cout << "process " << rank << " empty sync" << endl;
		ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
		ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
	}
	
	//Integrate IV curve
	//STILL UNCONFIRMED
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);
	if(rank == 0)
	{
		if(p.numEnergySteps >= 6)
		{
			for(int i=0; i<=p.vsdSteps; i++)
			{
				vsd = p.vsdMin + i*(p.vsdMax-p.vsdMin)/p.vsdSteps;
				corrsd = integrate(p, IVvalues, vsd);
				ierr = PetscPrintf(PETSC_COMM_SELF, "%.3f %.5e\n", vsd, corrsd);CHKERRQ(ierr);
				ierr = PetscFPrintf(PETSC_COMM_SELF, IVoutput, "%.3f %.5e\n", vsd, corrsd);CHKERRQ(ierr);
			} //endfor
		} //endif
		else
		{
			char message[] = "IV Integration not done. Must complete at least 6 energy steps for integration method to work.\n";
			ierr = PetscPrintf(PETSC_COMM_SELF, message);CHKERRQ(ierr);
			ierr = PetscFPrintf(PETSC_COMM_SELF, IVoutput, message);CHKERRQ(ierr);
		} //endelse
		ierr = PetscFClose(PETSC_COMM_WORLD, IVoutput);CHKERRQ(ierr);
	} //endif
	
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
	
	//Destroy PETSc Objects
    ierr = MatDestroy(identityMat);CHKERRQ(ierr);
    ierr = MatDestroy(inverseHamiltonian);CHKERRQ(ierr);
    ierr = MatDestroy(hamiltonian);CHKERRQ(ierr);
	ierr = KSPDestroy(ksp);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);

	//Destroy blocks (Is this necessary?  I believe class destructors are automatically executed when their scope ends.)
	for(int i=0; i<p.ndiag; i++)
		for(int j=0; j<3; j++)
			blockHamiltonian[i][j].~block();

	for(int i = 0; i < p.rows; i++)
	{
		free(*(atomCoordinates + i));
		free(*(disCoordinates + i));
	}
	free(atomCoordinates);
	free(disCoordinates);
	
	return 0;
} //end main
