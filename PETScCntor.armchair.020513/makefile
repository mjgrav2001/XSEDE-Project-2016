
CFLAGS 	        = 
FFLAGS	        = 
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = $PWD
EXAMPLESC       = cntor.cpp test.cpp ex1.c
EXAMPLESF       =
MANSEC          = Mat

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

cntor: cntor.o  chkopts
	-${CLINKER} -O0 -g -o cntor cntor.o ${PETSC_MAT_LIB}
	${RM} cntor.o

test: test.o  chkopts
	-${CLINKER} -O0 -g -o test test.o ${PETSC_MAT_LIB}
	${RM} test.o

ex1: ex1.o  chkopts
	-${CLINKER} -O0 -g -o ex1 ex1.o ${PETSC_MAT_LIB}
	${RM} ex1.o

runcntor: 
	-@${MPIEXEC}  -n 2 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps

rundebug: 
	-@${MPIEXEC}  -n 2 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -start_in_debugger -malloc_log

runarmchair:
	-@${MPIEXEC}  -n 2 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 100

runtest: 
	-@${MPIEXEC}  ./test

runex1: 
	-@${MPIEXEC}  ./ex1 

include ${PETSC_DIR}/conf/test
