/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
//#include <omp.h>

#include "array_2d.h"
#include "wave.h"
#include "array_2d_io.h"

int main(int argc, char** argv){

    if(argc != 6){
        printf("Incorrect number of parameters.  Correct usage:\n./wave_timing n Mx My alpha T\n");
        return 1;
    }
    MPI_Init(&argc, &argv);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Setup the timers.
    double start, stop;
    double simulation_time = 0.0;

    // Extract command line arguments.
    unsigned int n = atoi(argv[1]);
    unsigned int Mx = atoi(argv[2]);
    unsigned int My = atoi(argv[3]);

    float alpha = atof(argv[4]);
    int nt = atoi(argv[5]);

    // Perform the task.
//    start = omp_get_wtime();
    start = MPI_Wtime();
    standing_wave_simulation(nt, n, Mx, My, 1, alpha);

//    stop = omp_get_wtime();
    stop = MPI_Wtime();
    simulation_time = (double)(stop - start);

    printf("%d %d %d %f\n", size, n, nt, simulation_time);
    MPI_Finalize();
    return 0;
}
