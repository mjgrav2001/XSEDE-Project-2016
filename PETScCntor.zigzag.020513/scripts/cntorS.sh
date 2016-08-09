#!/bin/bash
#BSUB -J cntor
#BSUB -o stdout.txt
#BSUB -e stderr.txt
#BSUB -W 10
#BSUB -q small
#BSUB -n 1
#
./cntor > output.txt 2>&1
