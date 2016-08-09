#!/bin/csh
#BSUB -J cntor
#BSUB -o stdout.txt
#BSUB -e stderr.txt
#BSUB -a mpich2
#BSUB -W 24:00
#BSUB -q small
#BSUB -n 32
#
mpirun.lsf ./cntor lu
