#!/bin/bash
#BSUB -J cntor
#BSUB -o /nethome/mjack/Bench/PETSC/32/stdout.txt
#BSUB -e stderr.txt
#BSUB -a mpich2
#BSUB -W 6:00
#BSUB -q small
#BSUB -n 32
#
mpirun.lsf ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300