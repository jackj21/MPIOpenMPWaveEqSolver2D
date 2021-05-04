#!/bin/bash

#SBATCH --nodes=2
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

echo "wave images p=64"
echo " n Mx My" 
echo "wave error p=64"
echo 

mpirun -np 64 ./wave_images 1001 17 27 1.0 500
mpirun -np 64 ./wave_error 1001 17 27 1.0 500
