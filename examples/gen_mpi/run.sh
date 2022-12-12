#!/bin/bash 
#SBATCH -n 16
#SBATCH -t 1:00:00

mpprun ../../generator/build.mpi/bin/jgain_mpi.exe