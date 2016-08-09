#!/bin/bash
#BSUB -J cntor
#BSUB -o stdout.txt
#BSUB -e stderr.txt
#BSUB -a mpich2
#BSUB -W 10
#BSUB -q small
#BSUB -n 2
#
mpirun.lsf ./cntor -pc_type lu -pc_factor_mat_solver_package mumps