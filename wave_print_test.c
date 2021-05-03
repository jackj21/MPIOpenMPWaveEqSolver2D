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
#include <math.h>

#include "array_2d.h"
#include "wave.h"
#include "array_2d_io.h"

int main(int argc, char** argv){

    if(argc != 6){
        printf("Incorrect number of parameters.  Correct usage:\n./wave_print_test n Mx My alpha dt\n");
        return 1;
    }
    MPI_Init(&argc, &argv);
    
    // Extract command line arguments.
    unsigned int n = atoi(argv[1]);
    unsigned int Mx = atoi(argv[2]);
    unsigned int My = atoi(argv[3]);

    float alpha = atof(argv[4]);
    int nt = atoi(argv[5]);
	
	// Specificiation allows us to assume that nx = ny
	unsigned int nx = n;
	unsigned int ny = n;

	// Compute grid spacing
	float dx = 1.0 / (nx-1);
	float dy = dx;

	// Compute time step size
	float dt = alpha * dx / sqrt(2.0);

	int error = 0;

	Array2D_f u_prev;
	Array2D_f u_curr;
	Array2D_f u_next;
		
	allocate_Array2D_f(&u_prev, ny, nx, 0, MPI_COMM_WORLD);
	allocate_Array2D_f(&u_curr, ny, nx, 0, MPI_COMM_WORLD);
	allocate_Array2D_f(&u_next, ny, nx, 0, MPI_COMM_WORLD);

	initialize_Array2D_f(&u_prev);
	initialize_Array2D_f(&u_curr);
	initialize_Array2D_f(&u_next);
		
	// Evaluate initial conditions at t=-dt and t=0
	evaluate_standing_wave(&u_prev, Mx, My, dx, dy, -1.0*dt);
	evaluate_standing_wave(&u_curr, Mx, My, dx, dy, 0.0);

	char outfile[50];	
	
	int ny_local = u_curr.N_local;
	int nx_local = u_curr.N_local;
	int ny_padded = u_curr.N_padded;
	int nx_padded = u_curr.N_padded;
	int ny_global = u_curr.N_global;
	int nx_global = u_curr.N_global;
	int padding_ny = u_curr.padding;
	int padding_nx = u_curr.padding;
	float* arr = u_curr.data;
	MPI_Comm comm = u_curr.comm;


	memset(outfile, 50*sizeof(char), 0);
	sprintf(outfile, "wave_print_test_%05d.arr", 0);
	
	write_float_array_dist_cio(arr, ny_local, nx_local, ny_padded, nx_padded, ny_global, nx_global, padding_ny, padding_nx, comm, outfile);
	
	for (int k=1; k<nt; k++) {
		standing_wave_simulation_nsteps(&u_prev, &u_curr, &u_next, dt, dx, 1);

		memset(outfile, 50*sizeof(char), 0);
		sprintf(outfile, "wave_print_test_%05d.arr", k);
		write_float_array_dist_cio(arr, ny_local, nx_local, ny_padded, nx_padded, ny_global, nx_global, padding_ny, padding_nx, comm, outfile);
	}

	// Clean up memory and close MPI
	deallocate_Array2D_f(&u_prev);
	deallocate_Array2D_f(&u_curr);
	deallocate_Array2D_f(&u_next);

    MPI_Finalize();
    return 0;
}
