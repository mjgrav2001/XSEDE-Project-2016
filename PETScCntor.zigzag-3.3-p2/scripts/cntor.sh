#!/bin/bash
#MOAB -l nodes=2
#MOAB -o stdout.txt
#MOAB -e stderr.txt
#MOAB -j oe
#MOAB -l walltime=120:00
#MOAB -N cntor
cd $PBS_O_WORKDIR
mpirun ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -log_summary 2>&1 | output.txt