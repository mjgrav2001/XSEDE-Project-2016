#PBS -S /bin/bash
#PBS -q titan
#PBS -N cntor-zigzag-petsc-job
#PBS -j oe
#PBS -M mark.a.jack@gmail.com
#PBS -l walltime=1:00:00
#PBS -l nodes=1024
#PBS -A MAT101
#
source $MODULESHOME/init/bash
#module swap PrgEnv-pgi PrgEnv-intel
#module swap craype-interlagos craype-target-native
#
module load PrgEnv-intel
module load intel
module load cray-mpich
#
NNODES=1024
NMPITASKS=16
export OMP_NUM_THREADS=16
export MKL_MIC_ENABLE=1
export MIC_OMP_NUM_TRHREADS=96
#
PETSC_DIR=/ccs/home/n8d/petsc-3.5.0
#PETSC_DIR=$MEMBERWORK/mat101/petsc-3.5.0
export PETSC_DIR
#
#PETSC_ARCH=complex-cpp-mumps-titan
PETSC_ARCH=complex-cpp-optimized-titan
export PETSC_ARCH
#
#cd /ccs/home/n8d/ORNL-VFP-Project/cntor-petsc3.5/zigzag
#cd $PBS_O_WORKDIR
cd $MEMBERWORK/mat101/ORNL-VFP-Project/cntor-petsc3.5/zigzag/zigzag-10-0-30000/
#cd $PWD
#
aprun -n ${NNODES} -N 1 -d ${NMPITASKS} ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 arg1 2>&1 | tee /ccs/home/n8d/cntor-zigzag-N1024-pn16-$PBS_JOBID.log
#aprun -n ${NNODES} -N 1 -d ${NMPITASKS} ./cntor-petsc -pc_type lu arg1 2>&1 | tee /ccs/home/n8d/cntor-zigzag-petsc-N1024-pn1-30000-$PBS_JOBID.log