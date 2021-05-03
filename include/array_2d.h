/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#ifndef __ARRAY_2D_H__
#define __ARRAY_2D_H__

#include <mpi.h>

typedef struct Array2D_f_tag {

	//Local subset of data
    float* data;

    // ny X nx is number of rows X number of columns
    unsigned int ny;
    unsigned int nx;

	//Size of ghost padding
	int padding;

	//Size of (unpadded) local part of vector
	int N_local;

	//Size of padded local part of vector
	int N_padded;

	//Size of the global vector
	unsigned int N_global;

	//Starting idx of unpadded local vector in global coordinates
	int r0;

	//Communicator over which the global array is partitioned
	MPI_Comm comm;



} Array2D_f;

int allocate_Array2D_f(Array2D_f* arr, unsigned int ny, unsigned int nx, int padding, MPI_Comm comm);
int deallocate_Array2D_f(Array2D_f* arr);
int nullify_Array2D_f(Array2D_f* arr);
int initialize_Array2D_f(Array2D_f* arr);
int copy_Array2D_f(Array2D_f* src, Array2D_f* dest);
int write_Array2D_f(Array2D_f* u, char* filename);
int halo_exchange_Array2D(Array2D_f* arr);


#endif
