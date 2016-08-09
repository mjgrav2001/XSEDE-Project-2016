#!/bin/bash
#PBS -l nodes=8:ppn=8:procs8,pmem=2000M,walltime=16:00:00
#PBS -N cntor_armchair_job
#PBS -q batch128
#PBS -A mjack
#PBS -o stdout_armchair.log
#PBS -e error_armchair.log
#PBS -j oe
#PBS -V
 
# Enter Working Directory
cd $PBS_O_WORKDIR
export NP=`wc -l $PBS_NODEFILE | cut -d/ -f1`
export JOBID=`echo $PBS_JOBID | cut -d'.' -f1`
 
# Create log file for debugging
echo $PBS_O_WORKDIR
echo "Number of Processes: $NP" >> job_$JOBID.log
echo "JobID: $JOBID" >> job_$JOBID.log
echo "PBS NodeFile: $PBS_NODEFILE" >> job_$JOBID.log
cat $PBS_NODEFILE | uniq >> job_$JOBID.log
echo " " >> job_$JOBID.log
echo $SHELL >> job_$JOBID.log
echo " " >> job_$JOBID.log
 
# Run job
mpirun -np $NP ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 >> job_output.log