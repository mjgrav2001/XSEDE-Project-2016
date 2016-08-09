#!/bin/bash
#BSUB -J cntor
#BSUB -o stdout.txt
#BSUB -e stderr.txt
#BSUB -a mpich2
#BSUB -W 4:00
#BSUB -q small
#BSUB -n 64
#
#mpirun.lsf ./cntor -pc_type lu -pc_factor_mat_solver_package mumps
mpirun.lsf ./cntor
