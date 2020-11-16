//This program will work only when no. of processes = 10

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

float dot_product (void);

float a[2][8], b[8][10], c0, c1, d[2][10], parcel[8];
int my_rank;       /* rank of this process*/
int p;             /* number of processes*/

int main(int argc, char* argv[]) 
{
     MPI_Status status; /* return status*/
     
     int i, j, k;
     
     MPI_Init(&argc, &argv);//Start up MPI
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);//Find out the rank of this process
     MPI_Comm_size(MPI_COMM_WORLD, &p);//Find out number of processes
     
     if (my_rank == 0) {//initialize the matrices
          for (i = 0; i < 2; i++) {
               for (j = 0; j < 8; j++) {
                    a[i][j] = i+1.0;
               }
          }
          for (i = 0; i < 8; i++) {
               for (j = 0; j < 10; j++) {
                    b[i][j] = i*10.0 + j;
                    //b[i][j] = 1.0;
               }
          }
          /*
          //clear d
          for (i = 0; i < 2; i++) {
               for (j = 0; j < 10; j++) {
                    d[i][j] = 0.0;
               }
          }
          
          printf ("======Sequential===========\n");
          for (k = 0; k < 10; k++) {
               for (i = 0; i < 2; i++) {
                    for (j = 0; j < 8; j++) {
                         d[i][k] += a[i][j] * b[j][k];
                    }
               }
          }
          
          for (i = 0; i < 2; i++) {
               for (k = 0; k < 10; k++) {
                    printf ("%i%i:%f\t", i, k, d[i][k]);
               }
               printf ("\n"); 
          }
          */
     }

     //clear d again
     for (i = 0; i < 2; i++) {
          for (j = 0; j < 10; j++) {
               d[i][j] = 0.0;
          }
     }
     
     MPI_Bcast (a[0], 8, MPI_FLOAT, 0, MPI_COMM_WORLD);//broadcast matrix a to all nodes
     MPI_Bcast (a[1], 8, MPI_FLOAT, 0, MPI_COMM_WORLD);//broadcast matrix a to all nodes

     if (my_rank == 0) {//root does this
          for (k = 1; k < 10; k++) {//loop over all the nodes
               for (i = 0; i < 8; i++) {
                    //put each column of matrix b in a parcel
                    parcel[i] = b[i][k];
               }
               //send each column to node=k
               //int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
               MPI_Send(parcel, 8, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
          }
          //column 0 remains in node=0
          for (i = 0; i < 8; i++) {
               parcel[i] = b[i][0];
          }
     }
     else {//all nodes receive the parcel from root
          //int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
          MPI_Recv (parcel, 8, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
     }

     dot_product ();//compute the partial result in each node

     //int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
     MPI_Gather (&c0, 1, MPI_FLOAT, d[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
     MPI_Gather (&c1, 1, MPI_FLOAT, d[1], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

     if (my_rank == 0) {
          printf ("\n======Parallel===========\n");
          for (i = 0; i < 2; i++) {
               for (k = 0; k < 10; k++) {
                    printf ("%i%i:%.2f ", i, k, d[i][k]);
               }
               printf ("\n"); 
          }
     }

     MPI_Finalize();//shut down MPI
}

float dot_product (void)
{
     int i, j;
     
     c0 = c1 = 0.0;
     for (i = 0; i < 8; i++) {
          c0 += a[0][i] * parcel[i];
          c1 += a[1][i] * parcel[i];
     }
}