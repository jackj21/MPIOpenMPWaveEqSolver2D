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

echo "strong wave_timing"
echo "===================="
echo "p n nt time" > wave_time_strong.csv

	mpirun -np 1  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 2  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 4  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 8  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 16  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 32  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 64  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 128  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 256  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv
	mpirun -np 512  ./wave_timing 1001 17 27 1.0 25 >> wave_time_strong.csv

echo ""
echo ""
