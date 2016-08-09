#include "petscksp.h"
#include "petscmat.h" 
#include "classes.h"
#include "general.h"
#include "invert.h"
#include "grid_jobs.h"
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
	Parameters p(/*rows*/ 12000,				  
		     /*columns*/ 12000,
		     /*idim*/ 40,
		     /*energyMin*/ -0.018,
		     /*energyMax*/ 0.1,
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
		     /*groupSize*/ 4);
	
	//int *ranksInMyGroup;
	//ranksInMyGroup = new int[p.groupSize];
	//for(int l=0; l<p.groupSize; l++)
	//	ranksInMyGroup[l] = p.myGroupNum*p.groupSize + l;
	
	//double p.BField = (8./32.)*(0.0255664/2.)
	
	//int totalSteps = ((p.energyMin <= 0 && p.energyMax >= 0) ? p.numEnergySteps + 1 : p.numEnergySteps);
	
	int *isl, *isr, *isl180, *isr180;	
	isl = new int[p.leadSize];
	isr = new int[p.leadSize];
	isl180 = new int[p.leadSize];
	isr180 = new int[p.leadSize];
	
	isl[0] = 1;
	isl[1] = 2;	
	isl[2] = p.imax1;
	isl[3] = p.imax2;
	
	isr[0] = p.nrngd2;
	isr[1] = p.nrngd1;
	isr[2] = p.imaxn2;
	isr[3] = p.imaxn1;
	
	isl180[0] = 1;
	isl180[1] = 2;	
	isl180[2] = p.imax1;
	isl180[3] = p.imax2;
	
	isr180[0] = p.nrngd2;
	isr180[1] = p.nrngd1;
	isr180[2] = p.imaxn2;
	isr180[3] = p.imaxn1;
	
		
	//const int leadSize = 4;
	//int isl[leadSize] = {1, 2, p.imax1, p.imax2};
	//int isr[leadSize] = {p.imaxd21, p.imaxd22, p.imaxd21i, p.imaxd22i};
	//int isl180[leadSize] = {1, 2, p.imax1, p.imax2};
	//int isr180[leadSize] = {p.imaxd21, p.imaxd22, p.imaxd21i, p.imaxd22i};
	
	double atomCoordinates[p.rows][3];
	double disCoordinates[p.rows][3];
	
	double corrsd, vsd;
	double densityOfStates;
	double IVvalues[p.numEnergySteps];

	//Print run info
	if(p.worldSize % p.groupSize != 0 && p.numGroups > 1)
		PetscPrintf(PETSC_COMM_WORLD, "\n***WARNING - Non-optimal group configuration (last group is a partial group)***\n\n");
		
	if(p.myWorldRank == 0) {
		//cout << "rows " << p.rows << " idim " << p.idim << " ndiag " << p.ndiag << " numEnergySteps " << p.numEnergySteps << endl;
		ierr = PetscPrintf(PETSC_COMM_WORLD, "\nStart time: %s", ctime(&programStart));CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Total Processors: %d\n", p.worldSize);CHKERRQ(ierr);
		cout << "Groups: " << p.numGroups << endl;
		cout << "Group Size: " << p.groupSize << endl;
	}

	//Set up processor groups
	MPI_Comm PETSC_COMM_SUBGROUP;
	MPI_Comm_split(PETSC_COMM_WORLD, p.myGroupNum, p.myGroupRank, &PETSC_COMM_SUBGROUP);

	//Open output files
	FILE *outFile, *IVoutput;
	ierr = PetscFOpen(PETSC_COMM_WORLD,"DETEfEoutput.txt", "w", &outFile);CHKERRQ(ierr);
	ierr = PetscFOpen(PETSC_COMM_WORLD,"IVoutput.txt", "w", &IVoutput);CHKERRQ(ierr);
	
	if(p.myWorldRank == 0  && (outFile == NULL || IVoutput == NULL)) {
		cout << "Error: Output files not opened!\n";
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
	ReadInCoordinates(p.rows, atomCoordinates);	
	ReadInDisplaced(p.rows, disCoordinates);	

	//******************************************************************************************************
	//Setup PETSc Identity
	ierr = MatCreate(PETSC_COMM_SUBGROUP, &identityMat);CHKERRQ(ierr);
	ierr = MatSetType(identityMat, MATDENSE);CHKERRQ(ierr);	
	ierr = MatSetSizes(identityMat, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr);
	for(int i=0; i<p.rows; i++) {
		ierr = MatSetValue(identityMat, i, i, 1.0, INSERT_VALUES);CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Setup Hamiltonian
	ierr = MatCreate(PETSC_COMM_SUBGROUP, &hamiltonian);CHKERRQ(ierr);
	ierr = MatSetType(hamiltonian, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(hamiltonian, PETSC_DECIDE, PETSC_DECIDE, p.rows, p.columns);CHKERRQ(ierr);

	//Setup PETSc inverseMat
	ierr = MatDuplicate(identityMat, MAT_DO_NOT_COPY_VALUES, &inverseHamiltonian);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(inverseHamiltonian, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Setup PETSc Solver
	ierr = KSPCreate(PETSC_COMM_SUBGROUP, &ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, hamiltonian, hamiltonian, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
	ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	ierr = PCFactorSetZeroPivot(pc, eta1);CHKERRQ(ierr);
	//******************************************************************************************************

	//******************************************************************************************************
	//Create index vectors of contact sites between torus and metallic leads
	PetscInt numSubMatrices;
	static IS *indexSetRows, *indexSetColumns;
	
	int *vsl,*vsr;
	vsl = new int[p.leadSize];
	vsr = new int[p.leadSize];
	
	if(p.myGroupSize > 1) numSubMatrices = 2; //parallel
	else numSubMatrices = 1;

	
	for(int i=0; i < p.leadSize; i++) {
		vsl[i]=isl[i]-1;
		vsr[i]=isr[i]-1;
	}

	ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetRows);CHKERRQ(ierr);
	ierr = PetscMalloc(numSubMatrices*sizeof(IS **), &indexSetColumns);CHKERRQ(ierr);

	ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, indexSetRows);CHKERRQ(ierr);
	ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, indexSetColumns);CHKERRQ(ierr);
	
	if(numSubMatrices > 1) {
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsr, indexSetRows+1);CHKERRQ(ierr);
		ierr = ISCreateGeneral(PETSC_COMM_SUBGROUP, p.leadSize, vsl, indexSetColumns+1);CHKERRQ(ierr);
	}
	
	for(int i=0; i<numSubMatrices; i++) {
		ierr = ISSort(indexSetRows[i]);CHKERRQ(ierr);
		ierr = ISSort(indexSetColumns[i]);CHKERRQ(ierr);
	}

	delete[] vsl;
	delete[] vsr;
	
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
	for(p.energy = p.energyMin + p.myGroupNum*p.energyStep; p.energy < (p.energyMax + .9*p.energyStep); p.energy += p.energyStep * p.numGroups) {
		//Skipping energy = 0
		if(abs(p.energy) >= .9*p.energyStep)
		{
			ierr = wiregf(p, isl180, isr180, isl, isr, gl, gr, disCoordinates);CHKERRQ(ierr);
			ierr = InvertHamiltonian(PETSC_COMM_WORLD, p, hamiltonian, inverseHamiltonian, identityMat, pc, ksp, isl180, isr180, isl, isr, gl, gr, atomCoordinates, disCoordinates);CHKERRQ(ierr);
			ierr = RetrieveInformation(p, numSubMatrices, indexSetRows, indexSetColumns, inverseHamiltonian, isl, isr, densityOfStates, gl, gr, gaml, gamrr);CHKERRQ(ierr);
			ierr = CalculateTransport(PETSC_COMM_WORLD, p, densityOfStates, gl, gr, gaml, gamrr, IVvalues, outFile);CHKERRQ(ierr);
		}
		else {
			//Empty sync when skipping energy = 0 to ensure output order
			//cout << "rank " << p.myWorldRank << " empty sync at 0" << endl;
			ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
			ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
			ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
			ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
			//cout << "rank " << p.myWorldRank << " zero sync success" << endl;
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
	
	delete[] isl;
	delete[] isr;
	delete[] isl180;
	delete[] isr180;
	
	//delete[] ranksInMyGroup;

	//Destroy PETSc Objects
	for(int i=0; i<numSubMatrices; i++) {
		ierr = ISDestroy(indexSetRows[i]);CHKERRQ(ierr);
		ierr = ISDestroy(indexSetColumns[i]);CHKERRQ(ierr);
	}
	ierr = PetscFree(indexSetRows);CHKERRQ(ierr);
	ierr = PetscFree(indexSetColumns);CHKERRQ(ierr);

	//Destroy Matrices
	ierr = MatDestroy(identityMat);CHKERRQ(ierr);
	ierr = MatDestroy(hamiltonian);CHKERRQ(ierr);
	ierr = MatDestroy(inverseHamiltonian);CHKERRQ(ierr);
	ierr = KSPDestroy(ksp);CHKERRQ(ierr);
	ierr = PetscFClose(PETSC_COMM_WORLD, outFile);CHKERRQ(ierr);

	//Determine if extra step is needed for synchronized prints
	if((p.numEnergySteps+p.crossedZero) % p.numGroups != 0 && p.myGroupNum >= (p.numEnergySteps+p.crossedZero) % p.numGroups) {
		//cout << "rank " << p.myWorldRank << " empty sync at end" << endl;
		ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
		ierr = PetscSynchronizedFPrintf(PETSC_COMM_WORLD, outFile, "");CHKERRQ(ierr);
		ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);
		//cout << "rank " << p.myWorldRank << " empty sync success" << endl;
	}

	//ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);CHKERRQ(ierr);

	//Integrate IV curve
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);

	if(p.myWorldRank == 0) {
		if(p.numEnergySteps >= 6) {
			for(int i=0; i<=p.vsdSteps; i++) {
				vsd = p.vsdMin + i*(p.vsdMax-p.vsdMin)/p.vsdSteps;
				corrsd = integrate(p, IVvalues, vsd);
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
