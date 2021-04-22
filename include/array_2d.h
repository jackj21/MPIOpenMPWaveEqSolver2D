/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#ifndef __ARRAY_2D_H__
#define __ARRAY_2D_H__

typedef struct Array2D_f_tag {

    float* data;

    // ny X nx is number of rows X number of columns
    unsigned int ny;
    unsigned int nx;

} Array2D_f;

int allocate_Array2D_f(Array2D_f* arr, unsigned int ny, unsigned int nx);
int deallocate_Array2D_f(Array2D_f* arr);
int nullify_Array2D_f(Array2D_f* arr);
int initialize_Array2D_f(Array2D_f* arr);
int copy_Array2D_f(Array2D_f* src, Array2D_f* dest);
int write_Array2D_f(Array2D_f* u, char* filename);


#endif