#!/bin/bash
#MOAB -l nodes=1
#MOAB -j oe
#MOAB -l walltime=120:00
#MOAB -N test
cd $PBS_O_WORKDIR
mpiexec -n 2 ./test