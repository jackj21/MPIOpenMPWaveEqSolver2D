/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "array_2d.h"
#include "wave.h"

int main(int argc, char** argv){

    if(argc < 6){
        printf("Incorrect number of parameters.  Correct usage:\n./wave_images n Mx My alpha nt_plot0 nt_plot1 ...\n");
        return 1;
    }
	
	MPI_Init(&argc, &argv);

    // Setup the timers.
    double t_start, t_end;
    double t_total = 0.0;



    // Extract command line arguments.
    unsigned int n = atoi(argv[1]);
    unsigned int Mx = atoi(argv[2]);
    unsigned int My = atoi(argv[3]);

    float alpha = atof(argv[4]);

	MPI_Barrier(MPI_COMM_WORLD);

    // Perform the task.
    t_start = MPI_Wtime(); 

    // Specification allows us to assume that nx == ny
    unsigned int nx = n;
    unsigned int ny = n;

    // Compute grid spacing
    float dx = 1.0 / (nx - 1);
    float dy = dx;

    // Compute time step size
    float dt = alpha * dx / sqrt(2.0);

    int error = 0;

    Array2D_f u_prev;
    Array2D_f u_curr;
    Array2D_f u_next;
    Array2D_f u_true;

    // This was not required by the spec but it is an additional safety against using an
    // unallocated array.
    nullify_Array2D_f(&u_prev);
    nullify_Array2D_f(&u_curr);
    nullify_Array2D_f(&u_next);
    nullify_Array2D_f(&u_true);

    // Allocate the required arrays.
    allocate_Array2D_f(&u_prev, ny, nx, 1, MPI_COMM_WORLD);
    if (error) return 1;
    allocate_Array2D_f(&u_curr, ny, nx, 1, MPI_COMM_WORLD);
    if (error) return 1;
    allocate_Array2D_f(&u_next, ny, nx, 1, MPI_COMM_WORLD);
    if (error) return 1;
    allocate_Array2D_f(&u_true, ny, nx, 1, MPI_COMM_WORLD);
    if (error) return 1;

    // Initialize the required arrays.
    initialize_Array2D_f(&u_prev);
    initialize_Array2D_f(&u_curr);
    initialize_Array2D_f(&u_next);
    initialize_Array2D_f(&u_true);

    char outfile[50];

    // Loop over the remaining times and run the simulation to that time.
    // We have to reset every time because we don't enforce that the inputs are
    // sorted.  If we had them sorted, we could perform the simulation incrementally.
    for(int k=5; k<argc; k++){

        int nt = atoi(argv[k]);

        // Evaluate the initial conditions at t=-dt and t=0
        evaluate_standing_wave(&u_prev, Mx, My, dx, dy, -1.0*dt);    
        evaluate_standing_wave(&u_curr, Mx, My, dx, dy, 0.0);
        
        // Perform nt steps of the simulation
        standing_wave_simulation_nsteps(&u_prev, &u_curr, &u_next, dt, dx, nt);

        // Evaluate the true solution at that time point.  Note, this is not exactly
        // the same as the request final time
        evaluate_standing_wave(&u_true, Mx, My, dx, dy, nt*dt);
 
        // Save the simulated solutution
        memset(outfile, 50*sizeof(char), 0);
        sprintf(outfile, "computed_%05d.arr", nt);
        write_Array2D_f(&u_curr, outfile);

        // Save the true solution
        memset(outfile, 50*sizeof(char), 0);
        sprintf(outfile, "true_%05d.arr", nt);
        write_Array2D_f(&u_true, outfile);

        printf("Saving images from nt=%d, closest time is %f.\n", nt, nt*dt);
    }

    // Clean up the memory
    deallocate_Array2D_f(&u_prev);
    deallocate_Array2D_f(&u_curr);
    deallocate_Array2D_f(&u_next);
    deallocate_Array2D_f(&u_true);
    double stop = MPI_Wtime();
    double simulation_time = (double)(stop - t_start);
    printf("%f\n", simulation_time);
	MPI_Finalize();

    return 0;
}
