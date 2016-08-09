#!/bin/bash
#BSUB -J cntor
#BSUB -o /nethome/mjack/Bench/PETSC/32/zigzag32_out.txt
#BSUB -e stderr.txt
#BSUB -a mpich2
#BSUB -W 6:00
#BSUB -q small
#BSUB -n 32
#
mpirun.lsf ./cntor 