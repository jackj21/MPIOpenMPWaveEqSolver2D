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

    if(argc != 6){
        printf("Incorrect number of parameters.  Correct usage:\n./wave_animation n Mx My alpha T\n");
        return 1;
    }
	
	MPI_Init(&argc, &argv);


    // Extract command line arguments.
    unsigned int n = atoi(argv[1]);
    unsigned int Mx = atoi(argv[2]);
    unsigned int My = atoi(argv[3]);

    float alpha = atof(argv[4]);
    float T = atof(argv[5]);

    // Specification allows us to assume that nx == ny
    unsigned int nx = n;
    unsigned int ny = n;

    // Compute grid spacing
    float dx = 1.0 / (nx - 1);
    float dy = dx;

    // Compute time step size
    float dt = alpha * dx / sqrt(2.0);

    // Compute number of timesteps
    int nt = floor(T / dt);

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

    // Evaluate the initial conditions at t=-dt and t=0
    evaluate_standing_wave(&u_prev, Mx, My, dx, dy, -1.0*dt);    
    evaluate_standing_wave(&u_curr, Mx, My, dx, dy, 0.0);

    float* error_history = malloc((nt+1)*sizeof(float));

    for(int k=0; k<nt; ++k){

        evaluate_standing_wave(&u_true, Mx, My, dx, dy, k*dt);

        error_norm(&u_curr, &u_true, error_history+k);

        standing_wave_simulation_nsteps(&u_prev, &u_curr, &u_next, dt, dx, 1);
    }

    evaluate_standing_wave(&u_true, Mx, My, dx, dy, nt*dt);
    error_norm(&u_curr, &u_true, error_history+nt);

    char outfile[50];
    memset(outfile, 50*sizeof(char), 0);
    sprintf(outfile, "error_%08d.dat", n);

    FILE* f;
    f = fopen(outfile, "wb");
    int n_written;

    int to_write = nt+1;
    n_written = fwrite(&to_write, sizeof(unsigned int), 1, f);
    if (n_written != 1){
        fprintf(stderr, "Error writing ny to file %s: %d of 1 data written.\n", outfile, n_written);
        return 1;   
    }

    n_written = fwrite(error_history, sizeof(float), nt+1, f);
    if (n_written != nt+1){
        fprintf(stderr, "Error writing error history to file %s: %d of %d data written.\n", outfile, n_written, nt+1);
        return 1;   
    }

    fclose(f);
    free(error_history);

    // Clean up the memory
    deallocate_Array2D_f(&u_prev);
    deallocate_Array2D_f(&u_curr);
    deallocate_Array2D_f(&u_next);
    deallocate_Array2D_f(&u_true);

	MPI_Finalize();
}
