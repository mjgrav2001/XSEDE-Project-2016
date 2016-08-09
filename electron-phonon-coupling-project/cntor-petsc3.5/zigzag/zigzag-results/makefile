CFLAGS 		=
FFLAGS 		=
CPPFLAGS	= -Wall -O3 -openmp -I$(MKLROOT)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include/
#CPPFLAGS	= -Wall -O3 -openmp -DMKL_ILP64 -I$(MKLROOT)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include/boost 
FPPFLAGS        =
LDFLAGS         = 
LOCDIR          = $PWD
EXAMPLESC       = cntor.cpp test.cpp ex1.c
EXAMPLESF       =
MANSEC          = Mat


#Static linking:
#MKL_CUSTOM = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -lpthread -lm


#Dynamic linking:
MKL_CUSTOM = -L$(MKLROOT)/lib/intel64 -lmkl_intel64_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

all: cntor

#Normal, native compilation on login/service nodes
cntor-fortran: gfortran cntor-current-b0-0d32-blas-10-0-3600.f -openmp -o a.out ${MKL_CUSTOM}



#Normal, should tolerate most PETSc/MPI combinations
cntor: src/cntor.o chkopts
	-${CLINKER} -openmp -o cntor src/cntor.o ${PETSC_MAT_LIB} ${MKL_CUSTOM} ${BOOST}
	${RM} src/cntor.o


test: test.o  chkopts
	-${CLINKER}  -o test test.o ${PETSC_MAT_LIB}
	${RM} test.o

ex1: ex1.o  chkopts
	-${CLINKER}  -o ex1 ex1.o ${PETSC_MAT_LIB}
	${RM} ex1.o

runcntor: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps

rundebug: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -start_in_debugger

runarmchair: cntor
	-@${MPIEXEC}  -n 8 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 100

runtest: 
	-@${MPIEXEC}  ./test 

runex1: 
	-@${MPIEXEC}  ./ex1 

env-lib:
	@echo ${PETSC_MAT_LIB}

