#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=normal_q
#SBATCH --account=cmda3634_rjh
#SBATCH --exclusive
cd $SLURM_SUBMIT_DIR

#Load modules
module reset
module load gompi

make clean
make all

mpirun -np 1  ./wave_timing 2501 17 27 1.0 250
mpirun -np 2  ./wave_timing 5001 17 27 1.0 250
mpirun -np 4  ./wave_timing 10001 17 27 1.0 250
mpirun -np 8  ./wave_timing 20001 17 27 1.0 250
mpirun -np 16  ./wave_timing 40001 17 27 1.0 250
mpirun -np 32  ./wave_timing 80001 17 27 1.0 250
mpirun -np 64  ./wave_timing 160001 17 27 1.0 250
mpirun -np 128  ./wave_timing 320001 17 27 1.0 250
mpirun -np 256  ./wave_timing 640001 17 27 1.0 250
mpirun -np 512  ./wave_timing 1280001 17 27 1.0 250

