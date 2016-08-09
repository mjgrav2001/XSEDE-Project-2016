
CFLAGS 	        = 
FFLAGS	        = 
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = $PWD
EXAMPLESC       = cntor.cpp test.cpp test2.cpp test3.cpp ex1.c test2c.c
EXAMPLESF       =
MANSEC          = Mat

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

cntor: cntor.o  chkopts
	-${CLINKER}  -o cntor cntor.o ${PETSC_MAT_LIB}
	${RM} cntor.o

test: test.o  chkopts
	-${CLINKER}  -o test test.o ${PETSC_MAT_LIB}
	${RM} test.o

test2: test2.o  chkopts
	-${CLINKER}  -o test2 test2.o ${PETSC_MAT_LIB}
	${RM} test2.o

test2c: test2c.o  chkopts
	-${CLINKER}  -o test2c test2c.o ${PETSC_MAT_LIB}
	${RM} test2c.o

test3: test3.o  chkopts
	-${CLINKER}  -o test3 test3.o ${PETSC_MAT_LIB}
	${RM} test3.o

ex1: ex1.o  chkopts
	-${CLINKER}  -o ex1 ex1.o ${PETSC_MAT_LIB}
	${RM} ex1.o

runcntor: 
	-@${MPIEXEC}  ./cntor 

runtest: 
	-@${MPIEXEC}  ./test 

runex1: 
	-@${MPIEXEC}  ./ex1 

include ${PETSC_DIR}/conf/test
