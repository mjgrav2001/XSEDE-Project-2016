#!/bin/bash		  
#SBATCH -A TG-PHY140011         # Project name
#SBATCH -p normal	        # Queue name
#SBATCH -J cntor_armchair       # Job Name
#SBATCH -N 20                   # Requests 20 nodes total
#SBATCH -n 20		        # Requests 20/20=1 tasks per node
#SBATCH -c 16
#SBATCH -t 08:00:00		# Run time (hh:mm:ss) - 8.0 hours
#SBATCH --get-user-env
#SBATCH -o armchair.o%j.%N.%n   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -e armchair.e%j.%N.%n   # Name of the error file (eg. myMPI.eJobID)
#SBATCH --mail-type=all
#SBATCH --mail-user=mark.a.jack@gmail.com  # Address for email notification
#
#OMP_NUM_THREADS=1
#export OMP_NUM_THREADS
#export MKL_MIC_ENABLE=1
#export MIC_OMP_NUM_THREADS=96
#
#cd $WORK/ORNL_VFP_Project/cntor-petsc3.5/armchair/
cd $SLURM_SUBMIT_DIR
source ~/bin/tacc-modules.sh
#ibrun -n 20 -o 0 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 >> cntor_armchair_N1_n16.log # Run the MPI executable named "cntor"
ibrun ./cntor-petsc -n 24000 -pc_type lu >> cntor_armchair_petsc_10_10_24000_ph50.3067.log
