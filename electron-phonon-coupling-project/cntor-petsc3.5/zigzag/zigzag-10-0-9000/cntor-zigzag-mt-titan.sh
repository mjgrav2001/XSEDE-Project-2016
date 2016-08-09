#PBS -S /bin/bash
#PBS -V
#PBS -q titan
#PBS -N cntor-zigzag-job
#PBS -j oe
#PBS -M mark.a.jack@gmail.com
#PBS -l walltime=1:00:00
#PBS -l nodes=1
#PBS -A MAT101
#
NMPITASKS=1
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
#cd /ccs/home/n8d/ORNL-VFP-Project/cntor-petsc3.5
cd $PBS_O_WORKDIR
#
#aprun -n ${NMPITASKS} ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 arg1 2>&1 | tee /ccs/home/n8d/cntor-zigzag-N4-n16-$PBS_JOBID.log
aprun -n ${NMPITASKS} ./cntor-mt -pc_type lu arg1 2>&1 | tee /ccs/home/n8d/cntor-zigzag-mt-N4-n16-$PBS_JOBID.log
#aprun -n ${NMPITASKS} ./cntor-petsc -pc-type lu arg1 2>&1 | tee /ccs/home/n8d/cntor-zigzag-petsc-N4-n16-$PBS_JOBID.log
