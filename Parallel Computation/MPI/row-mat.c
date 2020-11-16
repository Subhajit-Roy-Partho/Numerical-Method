//This program will work only when the no. of processes = 5
//compile:
//mpicc row-mat.c

//run:
//mpiexec -n 5 ./a.out

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//function declaration
float dot_product (void);

//input matrices
float a[4], b[4][5];
//result matrix
float d[5];
//temporary variables used for computation
float c, parcel[4];

int main(int argc, char* argv[]) 
{
//* rank of this process
int my_rank;
//* total number of processes
int p;
//* return status
MPI_Status status;

int i, j;

     MPI_Init(&argc, &argv);//Start up MPI

     
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);//Find out the rank of this process
     
     MPI_Comm_size(MPI_COMM_WORLD, &p);//Find out the total number of processes

     //clear the result
     for (j = 0; j < 5; j++) {
          d[j] = 0.0;
     }
     
     if (my_rank == 0) {//initialize the vector and matrix
          for (i = 0; i < 4; i++) {
               a[i] = i + 1.0;
               for (j = 0; j < 5; j++) {
                    b[i][j] = i * 10 + j;
               }
          }
          /*
          //sequential multiplication
          printf ("======Sequential===========\n");
          for (j = 0; j < 5; j++) {
               for (i = 0; i < 4; i++) {
                    d[j] += a[i] * b[i][j];
               }
               printf ("%f\t", d[j]);
          }
          */
     }

     //again clear the result matrix
     for (j = 0; j < 5; j++) {
          d[j] = 0.0;
     }
     
     //broadcast the matrix a to all nodes
     MPI_Bcast (a, 4, MPI_FLOAT, 0, MPI_COMM_WORLD);//broadcast vector a to all nodes
     
     //this is done by the root process
     if (my_rank == 0) {
          for (j = 1; j < 5; j++) {
               //put each column into a temporary array
               for (i = 0; i < 4; i++) {
                    parcel[i] = b[i][j];
               }
               //send the column to node=j
               MPI_Send(parcel, 4, MPI_FLOAT, j, 0, MPI_COMM_WORLD);
          }
          //column 0 remains in node=0
          for (i = 0; i < 4; i++) {
               parcel[i] = b[i][0];
          }
     }
     else {//all other nodes receive the column from node=0
          MPI_Recv (parcel, 4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
     }
     
     //compute the partial result in all nodes, including root
     c = dot_product ();
     
     //collect the partial results from all nodes and put them into the result matrix in the root process
     MPI_Gather (&c, 1, MPI_FLOAT, d, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
     
     //print the result
     if (my_rank == 0) {
          printf ("\n======Parallel===========\n");
          for (j = 0; j < 5; j++) {
               printf ("%f\t", d[j]);
          }
          printf ("\n");
     }
     
     MPI_Finalize(); //shut down MPI
}

//this computes the partial resultin all nodes
float dot_product (void)
{
int i;
float c = 0.0;

     c = 0.0;
     for (i = 0; i < 4; i++) {
          c += a[i] * parcel[i];
     }
     
     return c;
}