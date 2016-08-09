#!/bin/bash
#PBS -A UT-BEACON-JACK
#PBS -l nodes=32,walltime=24:00:00
cd $PBS_O_WORKDIR
micmpiexec -n 512 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 arg1 2>&1 | tee cntor-armchair-n32-c4-512-$PBS_JOBID.log
