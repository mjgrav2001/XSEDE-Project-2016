#PBS -S /bin/bash
#PBS -V
#PBS -q titan
#PBS -N cntor-armchair-job
#PBS -j oe
#PBS -M mark.a.jack@gmail.com
#PBS -l walltime=1:00:00
#PBS -l nodes=4
#PBS -A MAT101
#
NCORES=16
export OMP_NUM_THREADS=1
export MKL_MIC_ENABLE=0
export MIC_OMP_NUM_TRHREADS=96
#
PETSC_DIR=/ccs/home/n8d/petsc-3.5.0
export PETSC_DIR
PETSC_ARCH=complex-ccp-mumps-titan
#PETSC_ARCH=complex-cpp-optimized-titan
export PETSC_ARCH
#
#cd /ccs/home/n8d/ORNL-VFP-Project-2014/cntor-petsc3.5
cd $PBS_O_WORKDIR
#
mpirun -v -np ${NCORES} -o 0 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 arg1 2>&1 | tee /ccs/home/n8d/cntor-armchair-N4-n16-$PBS_JOBID.log
#mpirun -v -np ${NCORES} -o 0 ./cntor-mt -pc_type lu arg1 2>&1 | tee /ccs/home/n8d/cntor-armchair-mt-N4-n16-$PBS_JOBID.log
#mpirun -v -np ${NCORES} -o 0 ./cntor-petsc -pc_type lu arg1 2>&1 | tee /ccs/home/n8d/cntor-armchair-petsc-N4-n16-$PBS_JOBID.log