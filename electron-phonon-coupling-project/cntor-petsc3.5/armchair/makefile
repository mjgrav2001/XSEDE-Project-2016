CFLAGS 	        = 
CXXFLAGS        =
FFLAGS	        = 
CPPFLAGS        = -Wall -O3 -openmp -I$(MKLROOT)/include -I$(TACC_BOOST_INC)
FPPFLAGS        =
LDFLAGS         = 
LOCDIR          = $PWD
EXAMPLESC       = cntor.cpp test.cpp ex1.c
EXAMPLESF       =
MANSEC          = Mat

BOOST = -L$(TACC_BOOST_LIB)
BOOST += -lboost_program_options

#Compile w/o MKL dense support
#CPPFLAGS += -DNOMKL

#Compile w/o boost libraries (must edit source for options)
#CPPFLAGS += -DNOBOOST

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

all: cntor

#Normal, should tolerate most PETSc/MPI combinations
cntor: src/cntor.o chkopts
	-${CLINKER} -o cntor src/cntor.o ${PETSC_MAT_LIB} ${BOOST}
	${RM} src/cntor.o

cntor-petsc: src/cntor.o chkopts
	-${CLINKER} -o cntor-petsc src/cntor.o ${PETSC_MAT_LIB} ${BOOST}
	${RM} src/cntor.o

runcntor: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps

rundebug: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -start_in_debugger

runarmchair: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 100

env:
	@echo "Libraries included:"
	@echo ${PETSC_MAT_LIB}
	@echo ${BOOST}
