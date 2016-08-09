#include "classes.h"
#include "general.h"
#include "grid_jobs.h"
#include "integration.h"
#include "invert.h"

#include <petscksp.h>
#include <petscmat.h>

#ifndef NOMKL
#include <omp.h>
#endif

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

	PetscInitialize(&argc, &args, NULL, NULL);
	PC pc;
	KSP ksp;
	PetscErrorCode ierr;
	Mat hamiltonian, inverseHamiltonian, identityMat;
	bool failFlag = false;

	//Numerical accuracies
	double eta1 = 1e-18;

	//Crtical values:
	//double bnatural = 0.004762;
	//Bfield = bnatural / 16.;

	//Setting all parameters in a struct for easier passing to functions and a central location to make changes
	//IMPORTANT: Subroutines will have their own *copy* of this class unless passed by reference.
	//			 If passed by reference, changes will be seen by all subsequent uses of the class.
	Parameters defaults(/*rows*/ 3600,
			    /*columns*/ 3600,
			    /*idim*/ 12,
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
			    /*a*/ 2.0,
			    /*groupSize*/ 1,
			    /*useDense*/ 0,
			    /*useZh*/ 0,
			    /*fileOptions*/ "w");
#ifndef NOBOOST
	Parameters p = applyCmdLine(defaults, argc, args);
#else
	Parameters p = defaults;
#endif

	int isl[p.leadSize];
	int isr[p.leadSize];
	int isl180[p.leadSize];
	int isr180[p.leadSize];
	
	isl[0] = 1;
	isl[1] = 2;	
	isl[2] = p.imax1;
	isl[3] = p.imax2;
	
	isr[0] = p.imaxd21;
	isr[1] = p.imaxd22;
	isr[2] = p.imaxd21i;
	isr[3] = p.imaxd22i;
	
	isl180[0] = 1;
	isl180[1] = 2;	
	isl180[2] = p.imax1;
	isl180[3] = p.imax2;
	
	isr180[0] = p.imaxd21;
	isr180[1] = p.imaxd22;	
	isr180[2] = p.imaxd21i;
	isr180[3] = p.imaxd22i;
	
	double atomCoordinates[p.rows][3];
	double disCoordinates[p.rows][3];
	
	double corrsd, vsd;
	double densityOfStates;
	//TODO eliminate this, since we obtain it from file
	double IVvalues[p.numEnergySteps];

	//Print run info
	if(p.worldSize % p.groupSize != 0 && p.numGroups > 1)
		PetscPrintf(PETSC_COMM_WORLD, "\n***WARNING - Non-optimal group configuration (last group is a partial group)***\n\n");
		
	if(p.myWorldRank == 0) {
		//cout << "rows " << p.rows << " idim " << p.idim << " ndiag " << p.ndiag << " numEnergySteps " << p.numEnergySteps << endl;
		cout << "\nStart time: " << ctime(&programStart) << endl;
		cout << "Total Processors: " << p.worldSize << endl;
		cout << "Groups: " << p.numGroups << endl;
		cout << "Processes per group: " << p.groupSize << endl;

#ifndef NOMKL
#pragma omp parallel
		if(omp_get_thread_num() == 0) cout << "Threads per process: " << omp_get_num_threads() << endl;
#endif
	}

	//Set up processor groups
	MPI_Comm PETSC_COMM_SUBGROUP;
	MPI_Comm_split(PETSC_COMM_WORLD, p.myGroupNum, p.myGroupRank, &PETSC_COMM_SUBGROUP);

	//Open output files
	FILE *outFile, *IVoutput;
	ierr = PetscFOpen(PETSC_COMM_WORLD,"DETEfEoutput.txt", p.fileOpts, &outFile);CHKERRQ(ierr);
	
	if(p.myWorldRank == 0  && outFile == NULL) {
		cout << "Error: Output file not opened!\n";
		failFlag = true;
	}
	if(failFlag) {
		ierr = PetscFinalize();CHKERRQ(ierr);
		return 0;
	}

	//Read in coordinates and setup system
	if(p.myWorldRank == 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Energy range: %.5f eV to %.5f eV\n", p.energyMin, p.energyMax);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, "BField = %f\n", p.BField);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\n Energy\t\tDensityofStates\t ScaledTransmissionFunction\t Re(fermif)\tTime\n\n");CHKERRQ(ierr);
	}
	//This takes about a second, since all procs access the same files
	ReadInCoordinates(p.rows, atomCoordinates);	
	ReadInDisplaced(p.rows, disCoordinates);	

	//******************************************************************************************************
	//Setup PETSc Identity for sparse solver, or not for dense
	if(!p.useDense) {
		ierr = MatCreate(PETSC_COMM_SUBGROUP, &identityMat);CHKERRQ(ierr);
		ierr = MatSetType(identityMat, MATDENSE);CHKERRQ(ierr);	//Might halve memory use?
		ierr = MatSetSizes(identityMat, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr); //+12000^2 * 8 = 1.1GB
		ierr = MatSetUp(identityMat);CHKERRQ(ierr);
		for(int i=0; i<p.rows; i++) {
			ierr = MatSetValue(identityMat, i, i, 1.0, INSERT_VALUES);CHKERRQ(ierr);
		}
		ierr = MatAssemblyBegin(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}

	//Setup Hamiltonian
	ierr = MatCreate(PETSC_COMM_SUBGROUP, &hamiltonian);CHKERRQ(ierr);
	if(p.useDense) {
		ierr = MatSetType(hamiltonian, MATSEQDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(hamiltonian, p.rows, p.columns, p.rows, p.columns);CHKERRQ(ierr);
	} else {
		ierr = MatSetType(hamiltonian, MATAIJ);CHKERRQ(ierr);
		ierr = MatSetSizes(hamiltonian, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr);
	}
	ierr = MatSetUp(hamiltonian);CHKERRQ(ierr); //TODO preallocate?

	if(p.useDense) {
		inverseHamiltonian = hamiltonian; //Dense done in-place
	} else {
		//Setup PETSc inverseHamiltonian
		ierr = MatDuplicate(identityMat, MAT_DO_NOT_COPY_VALUES, &inverseHamiltonian);CHKERRQ(ierr);
		ierr = MatSetUp(inverseHamiltonian);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		//Setup PETSc Solver
		ierr = KSPCreate(PETSC_COMM_SUBGROUP, &ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, hamiltonian, hamiltonian);CHKERRQ(ierr);
		//ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
		//ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
		ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
		ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
		ierr = PCFactorSetZeroPivot(pc, eta1);CHKERRQ(ierr);

		//Print solver type
		PCType type;
		ierr = PCGetType(pc, &type);CHKERRQ(ierr);
		if(p.myWorldRank == 0) cout << "PC type is: " << type << endl;
	}
	//******************************************************************************************************

	//******************************************************************************************************
	//Create index vectors of contact sites between torus and metallic leads
	PetscInt numSubMatrices;
	static IS *indexSetRows, *indexSetColumns;
	
	int vsl[p.leadSize];
	int vsr[p.leadSize];
	
	if(p.myGroupSize > 1) numSubMatrices = 2; //parallel
	else numSubMatrices = 1;

	
	for(int i=0; i < p.leadSize; i++) {
		vsl[i]=isl[i]-1;
		vsr[i]=isr[i]-1;
	}

	ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetRows);CHKERRQ(ierr);
	ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetColumns);CHKERRQ(ierr);

	ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, PETSC_USE_POINTER, indexSetRows);CHKERRQ(ierr);
	ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, PETSC_USE_POINTER, indexSetColumns);CHKERRQ(ierr);
	
	if(numSubMatrices > 1) {
	  ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, PETSC_USE_POINTER, indexSetRows+1);CHKERRQ(ierr);
	  ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, PETSC_USE_POINTER, indexSetColumns+1);CHKERRQ(ierr);
	}
	
	for(int i=0; i<numSubMatrices; i++) {
		ierr = ISSort(indexSetRows[i]);CHKERRQ(ierr);
		ierr = ISSort(indexSetColumns[i]);CHKERRQ(ierr);
	}
	
	//******************************************************************************************************
	//(step 0) declare arrays
	complex<double> **gl, **gr;
	complex<double> **gaml, **gamrr;
	
	//(step 1) allocate memory for array of elements of each column of coupling matrices
	gl = new complex<double>*[p.leadSize];
	gr = new complex<double>*[p.leadSize];
	gaml = new complex<double>*[p.leadSize];
	gamrr = new complex<double>*[p.leadSize];

	//(step 2) allocate memory for array of elements of each row of coupling matrices
	for( int i = 0; i < p.leadSize; i++) {
		*(gl + i) = new complex<double>[p.leadSize];
		*(gr + i) = new complex<double>[p.leadSize];
		*(gaml + i) = new complex<double>[p.leadSize];
		*(gamrr + i) = new complex<double>[p.leadSize];
	}

	//******************************************************************************************************
	//Loop over energy values
	//puts("Starting energy loop");
	for(p.energy = p.energyMin + p.myGroupNum*p.energyStep; p.energy < (p.energyMax + .9*p.energyStep); p.energy += p.energyStep * p.numGroups) {
		//Skipping energy = 0
		if(abs(p.energy) >= .9*p.energyStep)
		{
			ierr = wiregf(p, isl180, isr180, isl, isr, gl, gr, disCoordinates);CHKERRQ(ierr);
			ierr = InvertHamiltonian(PETSC_COMM_WORLD, p, hamiltonian, inverseHamiltonian, identityMat, pc, ksp, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates);CHKERRQ(ierr);
			//puts("Inverted Hamiltonian");
			ierr = RetrieveInformation(p, numSubMatrices, indexSetRows, indexSetColumns, inverseHamiltonian, isl, isr, densityOfStates, gl, gr, gaml, gamrr);CHKERRQ(ierr);
			ierr = CalculateTransport(PETSC_COMM_WORLD, p, densityOfStates, gl, gr, gaml, gamrr, IVvalues, outFile);CHKERRQ(ierr);
			//puts("Calculated transport");
		}
		else {
			//Empty sync when skipping energy = 0 to ensure output order
			ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
			ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, stdout);CHKERRQ(ierr);
			ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
			ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, outFile);CHKERRQ(ierr);
		}
	}
	//******************************************************************************************************

	//De-allocate memory
	for(int i=0; i<p.leadSize; i++) {
		delete[] *(gl+i);
		delete[] *(gr+i);
		delete[] *(gaml+i);	
		delete[] *(gamrr+i);			
	}
	delete[] gl;
	delete[] gr;
	delete[] gaml;	
	delete[] gamrr;

	//Destroy PETSc Objects
	for(int i=0; i<numSubMatrices; i++) {
	  ierr = ISDestroy(&(indexSetRows[i]));CHKERRQ(ierr);
	  ierr = ISDestroy(&(indexSetColumns[i]));CHKERRQ(ierr);
	}

	ierr = PetscFree(indexSetRows);CHKERRQ(ierr);
	ierr = PetscFree(indexSetColumns);CHKERRQ(ierr);

	//Destroy Matrices
	ierr = MatDestroy(&identityMat);CHKERRQ(ierr);
	ierr = MatDestroy(&hamiltonian);CHKERRQ(ierr);
	ierr = MatDestroy(&inverseHamiltonian);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = PetscFClose(PETSC_COMM_WORLD, outFile);CHKERRQ(ierr);

	//Determine if extra step is needed for synchronized prints
	if((p.numEnergySteps+p.crossedZero) % p.numGroups != 0 && p.myGroupNum >= (p.numEnergySteps+p.crossedZero) % p.numGroups) {
		//cout << "process " << p.myWorldRank << " empty sync" << endl;
		ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, stdout);CHKERRQ(ierr);
		ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, outFile);CHKERRQ(ierr);
	}

	//ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

	//Integrate IV curve
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);

	if(p.myWorldRank == 0) {

		ierr = PetscFOpen(PETSC_COMM_WORLD,"IVoutput.txt", "a", &IVoutput);CHKERRQ(ierr);

		double IVvals[p.numEnergySteps];
		ReadIVvalues(p,IVvals);

		if(p.numEnergySteps >= 6) {
			for(int i=0; i<=p.vsdSteps; i++) {
				vsd = p.vsdMin + i*(p.vsdMax-p.vsdMin)/p.vsdSteps;
				corrsd = integrate(p, IVvals, vsd);
				if(p.myGroupRank == 0) {
					ierr = PetscFPrintf(PETSC_COMM_WORLD, IVoutput, "%.3f %.5e\n", vsd, corrsd);CHKERRQ(ierr);
				}
			}
			ierr = PetscPrintf(PETSC_COMM_WORLD, "%.3f %.5e\n", vsd, corrsd);CHKERRQ(ierr);
		}
		else {
			char message[] = "IV Integration not done. Must complete at least 6 energy steps for integration method to work.\n";
			if(p.myGroupRank == 0) {
				ierr = PetscPrintf(PETSC_COMM_WORLD, message);CHKERRQ(ierr);
				ierr = PetscFPrintf(PETSC_COMM_WORLD, IVoutput, message);CHKERRQ(ierr);
			}
		}
		ierr = PetscFClose(PETSC_COMM_WORLD, IVoutput);CHKERRQ(ierr);
	}

	//Print run time
	if(p.myWorldRank == 0)
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

	ierr = MPI_Comm_free(&PETSC_COMM_SUBGROUP);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	
	return 0;

} //end main
