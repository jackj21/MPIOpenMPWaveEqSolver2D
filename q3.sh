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

echo "weak wave_timing"
echo "=================="
echo "p n nt time" > wave_time_weak.csv

for i in {1..3}
do
	mpirun -np 1  ./wave_timing 1001 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 2  ./wave_timing 2002 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 4  ./wave_timing 4004 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 8  ./wave_timing 8008 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 16  ./wave_timing 16016 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 32  ./wave_timing 32032 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 64  ./wave_timing 64064 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 128  ./wave_timing 128128 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 256  ./wave_timing 256256 17 27 1.0 250 >> wave_time_weak.csv
	mpirun -np 512  ./wave_timing 512512 17 27 1.0 250 >> wave_time_weak.csv
done
echo ""
echo ""
