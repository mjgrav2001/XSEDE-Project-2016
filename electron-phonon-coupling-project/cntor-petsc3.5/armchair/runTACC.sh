#!/bin/bash

export OMP_NUM_THREADS=16
export MKL_MIC_ENABLE=1
export MIC_OMP_NUM_THREADS=96
#cd $MAGMA_PATH/interface_mic/server
#./run.sh
#cd -
ibrun -n 16 -o 0 ./cntor -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 100
