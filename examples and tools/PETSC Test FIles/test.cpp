static char help[] = "Practice creating a matrix.\n";
#include <iostream>
#include "petscmat.h"
#undef __FUNCT__
#define __FUNCT__ "main"
using namespace std;

int main(int argc,char **args)
{
	PetscInitialize(&argc,&args,(char *)0,help);
	PetscErrorCode ierr;
	PetscInt size, rank;
	IS	isrow, iscol; 
	MatFactorInfo luinfo;
	Mat testMat, identityMat, inverseMat;
	PetscInt rows = 12, columns = 12, globalRowIDs[rows], globalColumnIDs[columns];
	PetscScalar identity[rows][columns], myMat[12][12] = {{10, 4,  6, 9,  9,  5, 1,  7, 2,  2, 10,  4},
														  { 0, 2,  5, 9,  3,  3, 2,  9, 6,  1,  9,  5},
														  { 5, 9,  5, 5,  9, 10, 8,  0, 9,  3,  6,  2},
												  		  { 9, 7, 10, 7,  6,  3, 4,  6, 2,  1,  1,  2},
														  { 3, 3,  2, 7,  3,  7, 7,  3, 9,  9,  1,  7},
														  { 7, 3,  4, 6,  4,  1, 7,  8, 8,  4,  6,  2},
														  { 7, 6,  9, 0, 10,  9, 2,  2, 7,  6,  7,  4},
														  { 6, 0,  9, 4,  8,  5, 3,  1, 8,  9, 10,  3},
														  { 7, 5,  8, 1,  0,  2, 7,  9, 1,  8,  1,  0},
														  { 9, 9,  9, 0,  8, 10, 1, 10, 7, 10,  8, 10},
													  	  { 0, 6,  2, 1,  1,  8, 4,  8, 6,  2,  3,  9},
														  { 4, 7,  2, 0,  7,  6, 0,  2, 5,  6,  7,  1}};

	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
	cout << "size: " << size << " rank: " << rank << endl;
	
	//Some Setup
	for(int i=0; i<rows; i++)
	{
		identity[i][i] = 1;
		globalRowIDs[i] = i;
		globalColumnIDs[i] = i;
	}
	
	//Setup testMat
	cout << "\ntestMat\n";
	ierr = MatCreate(PETSC_COMM_WORLD, &testMat);CHKERRQ(ierr); 
	ierr = MatSetSizes(testMat, PETSC_DECIDE, PETSC_DECIDE, rows, columns);CHKERRQ(ierr);
	ierr = MatSetType(testMat, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetValues(testMat, rows, globalRowIDs, columns, globalColumnIDs, &(myMat[0][0]), INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(testMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(testMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatView(testMat, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	//Setup identityMat. Note: duplicating a matrix is more efficient than making a new one.
	ierr = MatDuplicate(testMat, MAT_DO_NOT_COPY_VALUES, &identityMat);CHKERRQ(ierr);
	//ierr = MatSetType(identityMat, MATDENSE);CHKERRQ(ierr);
	ierr = MatSetValues(identityMat, rows, globalRowIDs, columns, globalColumnIDs, &(identity[0][0]), INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(identityMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	//Setup inverseMat
	ierr = MatDuplicate(testMat, MAT_DO_NOT_COPY_VALUES, &inverseMat);CHKERRQ(ierr);
	ierr = MatSetType(inverseMat, MATDENSE);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(inverseMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(inverseMat, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//Factor testMat and solve testMat * inverseMat = identityMat for inverseMat
    ierr = MatGetOrdering(testMat, MATORDERING_ND, &isrow, &iscol);CHKERRQ(ierr);
	ierr = MatFactorInfoInitialize(&luinfo);CHKERRQ(ierr);
	ierr = MatLUFactor(testMat, isrow, iscol, &luinfo);CHKERRQ(ierr);
	ierr = MatMatSolve(testMat, identityMat, inverseMat);CHKERRQ(ierr);
	
	//View Result. Note: the solution matrix, inverseMat, is the transpose of the true solution for some reason.
	cout << "\nfactoredMat\n";
	ierr = MatView(testMat, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	cout << "\ninverseMat\n";
	ierr = MatView(inverseMat, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	//Destroy Objects
    ierr = MatDestroy(testMat);CHKERRQ(ierr);
    ierr = MatDestroy(identityMat);CHKERRQ(ierr);
    ierr = MatDestroy(inverseMat);CHKERRQ(ierr);
	ierr = ISDestroy(isrow);CHKERRQ(ierr);
    ierr = ISDestroy(iscol);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
