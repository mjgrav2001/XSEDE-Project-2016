#!/bin/bash
#MOAB -l nodes=32
#MOAB -o Bench/PETSC/32/$PBS_JOBID.txt
#MOAB -m ae
#MOAB -e stderr.txt
#MOAB -j oe
#MOAB -l walltime=6:00:00
#MOAB -N cntor
cd $PBS_O_WORKDIR
mpirun ./cntor