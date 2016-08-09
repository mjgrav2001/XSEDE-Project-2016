#!/bin/bash
#MOAB -l nodes=16
#MOAB -o /panfs/storage.local/scs/home/mjack/$PBS_JOBID.armchair.txt
#MOAB -m ae
#MOAB -e stderr_$PBS_JOBID.txt
#MOAB -j oe
#MOAB -l walltime=16:00:00
#MOAB -N cntor
cd $PBS_O_WORKDIR
mpirun ./cntor