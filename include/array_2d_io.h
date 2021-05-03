#ifndef __ARRAY_IO_H__
#define __ARRAY_IO_H__

unsigned int subarray_size(unsigned int M, unsigned int P, int rank);

int write_float_array_dist_cio(float* arr,
                               int ny_local, int nx_local,
                               int ny_padded, int nx_padded,
                               int ny_global, int nx_global,
                               int padding_ny, int padding_nx,
                               MPI_Comm comm,
                               char* filename);

int write_float_array_dist_mpiio(float* arr,
                                 int ny_local, int nx_local,
                                 int ny_padded, int nx_padded,
                                 int ny_global, int nx_global,
                                 int padding_ny, int padding_nx,
                                 MPI_Comm comm,
                                 char* filename);

#endif
