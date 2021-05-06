#!/bin/bash

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=4
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

echo "all P wave_timing"
echo "==================="
echo "p n nt time" > wave_time_all_q5.csv
echo "p n nt time" > wave_time_all_q5_b.csv

for i in {1..3}
do
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 2501 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 5001 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 10001 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 25001 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 50001 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 75001 17 27 1.0 100 >> wave_time_all_q5.csv
	mpirun -np 512  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=2 ./wave_timing 100001 17 27 1.0 100 >> wave_time_all_q5.csv

	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 2501 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 5001 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 10001 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 25001 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 50001 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 75001 17 27 1.0 100 >> wave_time_all_q5_b.csv
	mpirun -np 256  --map-by ppr:32:node --bind-to L3cache -x OMP_NUM_THREADS=4 ./wave_timing 100001 17 27 1.0 100 >> wave_time_all_q5_b.csv



done

echo ""
echo ""
