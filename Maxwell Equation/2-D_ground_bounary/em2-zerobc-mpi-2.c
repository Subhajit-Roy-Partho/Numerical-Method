//em2-zerobc-mpi, version 2
//Parallel MPI 2D program with simple zero boundary conditions
//using (Ez, Hx, Hy)

//arrays via pointer allocation
//data saved via MPI-IO. Needs post processing.

//The mesh is spatially partitioned along x-direction
//The number of mesh points along the x-direction MUST BE divisible by the number of processors

//System dimensions are (maxi, maxj)
//Global array indices range from (0, 0) to (maxi-1, maxj-1)
//Internal points are from (1, 1) to (maxi-2, maxj-2)
//In each node mesh indices ranges from (node_mini, 0) to (node_maxi, maxj-1)

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

void exchange_boundary_fields (void);
long double **matrix (long rows, long cols);
void initialize (void);
void printinfo (char *text);

//sides of the box
int imax, jmax;
//limits within which E and H are computed in each node
int startE, startH, endE, endHx, endHy;

//variables required for mpi related stuff
//size of spatial partition in mesh units
int block_size, node_maxi, node_mini;
int my_rank, nproc;
MPI_Status status;
MPI_Offset offset, offset2;
MPI_File fh;//file handle

long double **ez, **hx, **hy;
long double pi, cee;
long double propfact;
int l, n, i, j, istart, jstart, nsteps;
long double dx, dt, t;
long double t0, spread, pulse;
FILE *fp;
char filename[20];

//storage for miscellaneous strings
char text[256];

int main (int argc, char** argv)
{
     MPI_Init(&argc, &argv);
     MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
     MPI_Comm_size( MPI_COMM_WORLD, &nproc );

     initialize ();
     
     nsteps = atoi (argv[1]);
     //delete all previous data files
     if (my_rank == 0)
          system ("rm Ez*");
     
     t = 0;
     
     for (n = 1; n <= nsteps ; n++) {
          pulse = exp (-0.5 * (pow ((t0 - t) / spread, 2.0)));
          
          //Start of the Main FDTD loop
          //Calculate the ez field
          for (j = 1; j < jmax-1; j++) {//1 to jmax-2
               for (i = startE; i <= endE; i++) {
                    ez[i][j] = ez[i][j] + propfact * (hy[i][j] - hy[i - 1][j] - hx[i][j] + hx[i][j - 1]) ;
               }
          }

          //Put a hard Gaussian pulse in the middle
          ez[istart][jstart] = pulse;
          
          //Calculate the hx field
          for (j = 0; j < jmax-1; j++) {
               for (i = startH; i <= endHx; i++) {
                    hx[i][j] = hx[i][j] + propfact * (ez[i][j] - ez[i][j + 1]) ;
               }
          }
          
          //Calculate the hy field
          for (j = 0; j < jmax; j++) {
               for (i = startH; i <= endHy; i++) {
                    hy[i][j] = hy[i][j] + propfact * (ez[i + 1][j] - ez[i][j]) ;
               }
          }
          
          //End of the main FDTD loop
          //synchronize all nodes
          MPI_Barrier(MPI_COMM_WORLD);
          
          //MPI-IO, all processes write different portions of the same file
          //NOTE: the file so created stores the data in binary, not text
          //filename is created at runtime from timestep
          sprintf (filename, "Ez.%i", n);
          MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
          
          //this will not work! ez is not a true array. See version 3
          //MPI_File_write_at(fh, offset, ez[node_mini], block_size*jmax, MPI_LONG_DOUBLE, &status);
          
          //but this will work, since each ez[i] points to a contiguous array
          j = 0;
          for (i = node_mini; i <= node_maxi; i++) {
               offset2 = j * jmax * sizeof (long double);
               MPI_File_write_at(fh, offset+offset2, ez[i], jmax, MPI_LONG_DOUBLE, &status);
               j++;
          }
          
          //close data file for this timestep 
          MPI_File_close(&fh);
          
          exchange_boundary_fields ();
          t = t + 1;
          
     }//for n - timestep ends
     
          
     MPI_Finalize();
     return 0;
     
}

//end main

void exchange_boundary_fields (void)
{
int i, j, k, l;
MPI_Status status;
     
     //ez sent from process l+1 and received in process l
     for (l = 0; l < nproc - 1; l++) {
          if (my_rank == l) {
               i = MPI_Recv (&(ez[node_maxi+1][0]), jmax, MPI_LONG_DOUBLE, l+1, 0, MPI_COMM_WORLD, &status);
               
          }
          if (my_rank == l+1) {
               MPI_Send (&(ez[node_mini][0]), jmax, MPI_LONG_DOUBLE, l, 0, MPI_COMM_WORLD);
               
          }
     }
     
     //hy transferred from process l-1 to l
     for (l = 1; l < nproc; l++) {
          if (my_rank == l) {
               MPI_Recv (&(hy[node_mini-1][0]), jmax, MPI_LONG_DOUBLE, l-1, 0, MPI_COMM_WORLD, &status);
               
          }
          if (my_rank == l-1) {//left process sends
               MPI_Send (&(hy[node_maxi][0]), jmax, MPI_LONG_DOUBLE, l, 0, MPI_COMM_WORLD);
               
          }
     }    
}

//end exchange_boundary_fields

void initialize (void)
{
//mesh is spatially partitioned along x-direction
//imax MUST BE divisible by nproc
     
int i, j;
long double p;
     
     pi = 3.14159;
     cee = 3.0e8;
     imax = 60;
     jmax = 30;
     
     //source is at the centre
     istart = imax / 2;
     jstart = jmax / 2;
     
     dx = 0.01; //Cell size
     dt = dx / (2.0 * cee); //Timestep
     
     //parameters for the gaussian pulse
     t0 = 80.0;
     spread = 12.0;
     
     //used for advancing Ez, Hx and Hy inside the system
     propfact = cee * dt / dx;
     
     //this is the size of the partition in each node along x-direction
     //it is the no. of mesh points along x in each node
     block_size = imax / nproc;
     
     //use rank of this node to find the the minimum and maximum of the global index along x-direction
     //in each node mesh index ranges from (node_mini, 0) to (node_maxi, maxj-1)
     node_mini = block_size * my_rank;
     node_maxi = node_mini + block_size - 1;
     
     //general case for all nodes
     startE = startH = node_mini;
     endE = endHx = endHy = node_maxi;
     
     if (my_rank == 0)//first node
          startE = 1;
     if (my_rank == nproc - 1) {//last node
          endE = endHy = node_maxi-1;
     }
     
     sprintf (text, "Block size %i Mini %i, Maxi %i", block_size, node_mini, node_maxi);
     printinfo (text);
     
     //offset to separate file write operation from different nodes
     offset = my_rank * block_size * jmax * sizeof (long double);
     
     ez = matrix (imax, jmax);
     hx = matrix (imax, jmax);
     hy = matrix (imax, jmax);
     
     for (i = node_mini; i <= node_maxi; i++) {
          for (j = 0; j < jmax; j++){
               ez[i][j] = 0.0;
               hx[i][j] = 0.0;
               hy[i][j] = 0.0;
          }
     }
     //synchronize all nodes
     MPI_Barrier(MPI_COMM_WORLD);

     //sprintf (text, "Init done.\n");
     //printinfo (text);
     
}

//end initialize

long double **matrix (long rows, long cols)
{
     int i;
     long double **p;
     
     p = (long double **) calloc (rows, sizeof(long double *));
     for (i = 0; i < rows; i=i+1) {
          p[i] = calloc(cols, sizeof (long double));
     }
     
     return p;
     
}

//end matrix

void printinfo (char *text)
{
     int i;
     MPI_Status status;
     
     if (my_rank == 0) {//root
          printf ("0: %s\n", text);//print own message
          for (i = 1; i < nproc; i++) {//receive messages from other nodes and print them
               MPI_Recv (text, 256, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
               printf ("%i: %s\n", i, text);
          }
     }
     else {//other nodes send to root
          MPI_Send (text, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
     }
     
     fflush (NULL);
}

//end printinfo

//end file
