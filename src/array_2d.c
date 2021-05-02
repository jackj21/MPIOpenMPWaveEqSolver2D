/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <mpi.h>
#include "array_2d.h"



/*
* Allocate an Array2D_f struct.
*
* Arguments:
*   arr: The ny x nx array struct to allocate into
*   m: number of rows in the array
*   n: number of columns in the array
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int allocate_Array2D_f(Array2D_f* arr, unsigned int m, unsigned int n, int padding, MPI_Comm comm) {

	int rank, size;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	int N_local, r0;
	//compute subvector size to find N_local and r0
	//look at "flag" in lab 8
	
	//Assign padding and size of local vector
	arr->padding = padding;
	arr->N_local = N_local;
	arr->N_padded = N_local + 2*padding;

	//Assign x and y dimensions
    arr->ny = m;
    arr->nx = n;
	//Assign size of global vector
    arr->N_global = arr->ny * arr->nx;

	//Assign the communicator (might need to change if using I/O)
	arr->comm = comm;
	
	//May not need the (float*)
	arr->data = (float*)malloc((arr->N_padded)*sizeof(float));

    if (arr->data == NULL){
        fprintf(stderr, "Error allocating 2D int array.\n");
        return 1;
    }

	//will return "flag" once compute_subarray_function completed
    return 0;
}


/*
* Deallocate an Array2D_f struct.
*
* Arguments:
*   arr: The ny x nx array struct to deallocate
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int deallocate_Array2D_f(Array2D_f* arr){

    arr->ny = 0;
    arr->nx = 0;
	arr->N_global = 0;

	arr->padding = 0;
	arr->N_local = 0;
	arr->N_padding = 0;

	arr->r0 = -1;

	arr->comm = MPI_COMM_NULL;

    free(arr->data);

    arr->data = NULL;

    return 0;
}


/*
* Initialize an Array2D_f struct to 0.
*
* Arguments:
*   arr: The ny x nx array struct to initialize
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int initialize_Array2D_f(Array2D_f* arr){

    memset(arr->data, 0, arr->N_padded*sizeof(float));

    return 0;
}


/*
* Nullify an Array2D_f struct.
*
* Nullification means pointers are set to NULL and dimensions are set to 0.
*
* Arguments:
*   arr: The ny x nx array struct to deallocate
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int nullify_Array2D_f(Array2D_f* arr){

    arr->ny = 0;
    arr->nx = 0;
	arr->N_global = 0;

	arr->N_padded = 0;
	arr->N_local = 0;
	arr->r0 = -1;
	arr->padding = 0;
	arr->comm = MPI_COMM_NULL;

    arr->data = NULL;

    return 0;
}


/*
* Copy an Array2D_f struct to another Array2D_f struct.
*
* Arguments:
*   src: The array to copy from
*   dest: The array to copy to
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int copy_Array2D_f(Array2D_f* src, Array2D_f* dest){

    if ((src->ny != dest->ny) ||
        (src->nx != dest->nx)){
        fprintf(stderr, "Size mismatch in int array copy (src=%dX%d dest=%d,%d).\n",
                src->ny, src->nx, dest->ny, dest->nx);
        return 1;
    }

    if (src->data == NULL){
        fprintf(stderr, "Source array NULL in int array copy.\n");
        return 2;
    }

    if (dest->data == NULL){
        fprintf(stderr, "Destination array NULL in int array copy.\n");
        return 3;
    }

    unsigned int n_data = src->ny * src->nx;

    memcpy(dest->data, src->data, n_data*sizeof(float));

    return 0;
}


/*
* Writes an Array2D_f struct to disk.
*
* Arguments:
*   u: The array to write
*   filename: The filename to write to
*
* Returns:
*   error code, 0 for success, 1 for failure
*/

int write_Array2D_f(Array2D_f* u, char* filename){

    FILE* f;
    unsigned char val;
    int idx;

    // Open the file
    f = fopen(filename, "wb");
    if (f == NULL) {
        fprintf(stderr, "Error opening file %s for write.\n", filename);
        return 1;
    }

    int n_written;

    // Write the header
    n_written = fwrite(&u->ny, sizeof(unsigned int), 1, f);
    if (n_written != 1){
        fprintf(stderr, "Error writing ny to file %s: %d of 1 data written.\n", filename, n_written);
        return 1;   
    }

    n_written = fwrite(&u->nx, sizeof(unsigned int), 1, f);
    if (n_written != 1){
        fprintf(stderr, "Error writing nx to file %s: %d of 1 data written.\n", filename, n_written);
        return 1;   
    }

    // Write the data

    unsigned int n_data = u->nx * u->ny;
    n_written = fwrite(u->data, sizeof(float), n_data, f);
    if (n_written != n_data){
        fprintf(stderr, "Error writing data to file %s: %d of %d data written.\n", filename, n_written, n_data);
        return 1;   
    }

    fclose(f);

    return 0;
}
