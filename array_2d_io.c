#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

unsigned int subarray_size(unsigned int M, unsigned int P, unsigned int rank){

//    <student>

    printf("You should remove this line!  If you see it in your output you didn't finish this function first!\n")

//    </student>

}

int write_float_array_dist_cio(float* arr,
                               int ny_local, int nx_local,
                               int ny_padded, int nx_padded,
                               int ny_global, int nx_global,
                               int padding_ny, int padding_nx,
                               MPI_Comm comm,
                               char* filename){

    // File handle
    FILE* f;

    // Total data to write, measured data written
    int n_data, n_written;

    // header size in bytes; two unsigned ints
    int header_size = 8;

    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Only rank 0 should manage file creation and write the header
    if (rank == 0){
        // Open the file
        f = fopen(filename, "wb");
        if (f == NULL) {
            fprintf(stderr, "Error opening file %s for write.\n", filename);
            return 1;
        }

        // Write the header
        n_written = fwrite(&ny_global, sizeof(unsigned int), 1, f);
        if (n_written != 1){
            fprintf(stderr, "Error writing ny to file %s: %d of 1 data written.\n", filename, n_written);
            return 1;   
        }

        n_written = fwrite(&nx_global, sizeof(unsigned int), 1, f);
        if (n_written != 1){
            fprintf(stderr, "Error writing nx to file %s: %d of 1 data written.\n", filename, n_written);
            return 1;   
        }
    }

    // The total data this rank will be writing
    n_data = ny_local*nx_local;

    // This is a little inefficient, but for simplicity we will create a local,
    // unpadded arrays so we can do a single communication.  This is safe because
    // by specification, Rank 0 has the largest possible array.
    float* write_data = (float*) malloc(n_data*sizeof(float));
    if(write_data == NULL){
        fprintf(stderr, "Error allocating temporary array data.\n");
        return 1;
    }
    memset(write_data, 0, n_data*sizeof(float));

    // Loop over the local, non-padded part of the data and fill
    // with the data
    for(int i=0;i<ny_local;i++){

        // Copy the ith row into the ith row of write_data.
        // The ith row starts padding_ny rows into the provided array
        // and the data we care about starts padding_nx columns into
        // that row.
        // The length of the row is nx_local.
        memcpy(write_data + i*nx_local, 
               arr + padding_ny*nx_padded + i*nx_padded + padding_nx,
               nx_local*sizeof(float));
    }

    if (rank == 0){

        // Rank 0 can directly write its own data
        n_written = fwrite(write_data, sizeof(float), n_data, f);
        if(n_written != n_data){
            fprintf(stderr, "Error writing file %s:expected %d bytes wrote %d bytes.\n", filename, n_data, n_written);
            return 1;
        }

        // Rank 0 receives chunks of data from all other processors and writes
        // it to disk
        int ny_p, nx_p;
        for(int p=1; p<size; p++){

            // Determine the remote ranks array properties (local rows and
            // columns, amount of data, etc.)
            ny_p = subarray_size(ny_global, size, p);
            nx_p = subarray_size(nx_global, 1, p);
            n_data = ny_p*nx_p;

            memset(write_data, 0, n_data*sizeof(float));

            // Rank 0 obtains part of the array from every other rank
            MPI_Status status;
            MPI_Recv(write_data, n_data, MPI_FLOAT, p,  0,  comm, &status);

            n_written = fwrite(write_data, sizeof(float), n_data, f);
            if(n_written != n_data){
                fprintf(stderr, "Error writing file %s:expected %d bytes wrote %d bytes.\n", filename, n_data, n_written);
                return 1;
            }
        }

        fclose(f);
    }
    else{
        // All other ranks just send their respective chunk of the array
        // to rank 0
        MPI_Send(write_data, n_data, MPI_FLOAT, 0,  0,  comm);
    }

    free(write_data);

    return 0;

}

int write_float_array_dist_mpiio(float* arr,
                                 int ny_local, int nx_local,
                                 int ny_padded, int nx_padded,
                                 int ny_global, int nx_global,
                                 int padding_ny, int padding_nx,
                                 MPI_Comm comm,
                                 char* filename){

    // File handle
    FILE* f;

    // Total data to write, measured data written
    int n_data, n_written;

    // header size in bytes; two unsigned ints
    int header_size = 8;

    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Only rank 0 should manage file creation and write the header
    if (rank == 0){
        // Open the file
        f = fopen(filename, "wb");
        if (f == NULL) {
            fprintf(stderr, "Error opening file %s for write.\n", filename);
            return 1;
        }

        // Write the header
        n_written = fwrite(&ny_global, sizeof(unsigned int), 1, f);
        if (n_written != 1){
            fprintf(stderr, "Error writing ny to file %s: %d of 1 data written.\n", filename, n_written);
            return 1;   
        }

        n_written = fwrite(&nx_global, sizeof(unsigned int), 1, f);
        if (n_written != 1){
            fprintf(stderr, "Error writing nx to file %s: %d of 1 data written.\n", filename, n_written);
            return 1;   
        }
        fclose(f);
    }

    // The total data this rank will be writing
    n_data = ny_local*nx_local;

    // This is a little inefficient, but for simplicity we will create a local,
    // unpadded arrays so we can do a single communication.  This is safe because
    // by specification, Rank 0 has the largest possible array.
    float* write_data = (float*) malloc(n_data*sizeof(float));
    if(write_data == NULL){
        fprintf(stderr, "Error allocating temporary array data.\n");
        return 1;
    }
    memset(write_data, 0, n_data*sizeof(float));

    // Loop over the local, non-padded part of the data and fill
    // with the data
    for(int i=0;i<ny_local;i++){

        // Copy the ith row into the ith row of write_data.
        // The ith row starts padding_ny rows into the provided array
        // and the data we care about starts padding_nx columns into
        // that row.
        // The length of the row is nx_local.
        memcpy(write_data + i*nx_local, 
               arr + padding_ny*nx_padded + i*nx_padded + padding_nx,
               nx_local*sizeof(float));
    }

    MPI_File out_file;
    MPI_Status status;
    MPI_Offset offset;

    // Perform the actual write using MPI I/O
    // <student>

    // </student>

    free(write_data);

    return 0;

}