#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define ARRAYSIZE 100

int main(int argc, char* argv[]) 
{
MPI_Status status;
int i, psize, my_rank, err;
int blocksize;

float *x, *y;
     
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
     MPI_Comm_size(MPI_COMM_WORLD, &psize);
     
     if (my_rank == 0) {
          x = (float *) calloc (ARRAYSIZE, sizeof (float));
          
          for (i = 0; i < ARRAYSIZE; i++) {
               x[i] = i * i;
               printf ("x[%i] = %f \n", i, x[i]);
          }
          printf ("\n");
     }

     blocksize = ARRAYSIZE/psize;
     y = (float *) calloc (blocksize, sizeof (float));
     
     err = MPI_Scatter(x, blocksize, MPI_FLOAT, y, blocksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
     
     for (i = 0; i < blocksize; i++) {
          y[i] = sqrt(y[i]);
     }
     if (my_rank == 1) {
          for (i = 0; i < blocksize; i++) {
               printf ("y[%i] = %f\n", i, y[i]);
          }
     }
     
     err = MPI_Gather (y, blocksize, MPI_FLOAT, x, blocksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
     if (my_rank == 0) {
          for (i = 0; i < ARRAYSIZE; i++) {
               //printf ("x[%i] = %f\n", i, x[i]);
          }
     }
     
     MPI_Finalize();//shut down MPI
}