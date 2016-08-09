#!/bin/bash
#PBS -A cnms14q
#PBS -N cntor-armchair-job
#PBS -j oe
#PBS -M mark.a.jack@gmail.com
#PBS -l walltime=2:00:00,nodes=8:ppn=8
NCORES=64
EXEC=./cntor
cd $PBS_O_WORKDIR
mpirun -v --mca m_pinned 1 --mca mpool_base_use_mem_hooks 1 --mca bt1 openib,self -np ${NCORES} ${EXEC} -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 arg1 2>&1 | tee cntor-armchair-n8-c8-64-$PBS_JOBID.log