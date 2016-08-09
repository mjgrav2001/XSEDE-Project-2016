#!/bin/bash		  
#SBATCH -A TG-PHY140011         # Project name
#SBATCH -p normal	        # Queue name
#SBATCH -J cntor_zigzag         # Job Name
#SBATCH -N 1                    # Requests 1 nodes total
#SBATCH -n 16		        # Requests 16/1=16 tasks per node
#SBATCH -t 01:00:00		# Run time (hh:mm:ss) - 1.0 hours
#SBATCH -o zigzag.o%j.%N.%n     # Name of the output file (eg. myMPI.oJobID)
#SBATCH -e zigzag.e%j.%N.%n     # Name of the error file (eg. myMPI.eJobID)
#SBATCH --mail-user=mark.a.jack@gmail.com  # Address for email notification
OMP_NUM_THREADS=16
export OMP_NUM_THREADS
export MKL_MIC_ENABLE=1
export MIC_OMP_NUM_THREADS=96
#cd $WORK/ORNL_VFP_Project/PETScCntor.zigzag.061214/
cd $SLURM_SUBMIT_DIR
ibrun -n 16 -o 0 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 300 >> cntor_zigzag_N1_n16.log # Run the MPI executable named "cntor"
#ibrun -n 16 -o 0 ./cntor-mkl >> cntor_mkl_zigzag_N1_n16.log
#ibrun -n 16 -o 0 ./cntor-petsc >> cntor_petsc_zigzag_N1_n16.log
#ibrun -n 16 -o 0 ./cntor-mt >> cntor_mt_zigzag_N1_n16.log
#ibrun -n 16 -o 0 ./cntor-mic >> cntor_mic_zigzag_N1_n16.log
