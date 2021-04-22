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

int allocate_Array2D_f(Array2D_f* arr, unsigned int m, unsigned int n) {

    arr->ny = m;
    arr->nx = n;

    unsigned int n_data = m*n;

    arr->data = malloc(n_data*sizeof(float));

    if (arr->data == NULL){
        fprintf(stderr, "Error allocating 2D int array.\n");
        return 1;
    }

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

    unsigned int n_data = arr->ny * arr->nx;

    memset(arr->data, 0, n_data*sizeof(float));

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
