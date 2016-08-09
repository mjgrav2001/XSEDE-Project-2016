#!/bin/bash
#MOAB -l nodes=32
#MOAB -o Bench/PETSC/32/$PBS_JOBID.txt
#MOAB -m ae
#MOAB -e stderr.txt
#MOAB -j oe
#MOAB -l walltime=16:00:00
#MOAB -N cntor
cd $PBS_O_WORKDIR
mpirun ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300