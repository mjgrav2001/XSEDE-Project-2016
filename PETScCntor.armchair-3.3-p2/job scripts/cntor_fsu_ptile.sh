#!/bin/bash
#MOAB -l nodes=64
#MOAB -o stdout.txt
#MOAB -e stderr.txt
#MOAB -j oe
#MOAB -l walltime=120:00
#MOAB -N cntor
#MOAB -R "span[ptile=8]"
cd $PBS_O_WORKDIR
mpirun ./cntor lu