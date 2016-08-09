#!/bin/bash
#BSUB -J cntor
#BSUB -o /nethome/mjack/zigzag64_10_0_12000_out.txt
#BSUB -e stderr_zigzag.txt
#BSUB -a mpich2
#BSUB -W 16:00
#BSUB -q medium
#BSUB -n 64
#
mpirun.lsf ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 1000