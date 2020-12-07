//em2-zerobc-mpi, version 1
//Parallel MPI 2D program with simple zero boundary conditions
//using (Ez, Hx, Hy)

//arrays via pointer allocation
//data saved separately by each process

//The mesh is spatially partitioned along x-direction
//The number of mesh points along the x-direction MUST BE divisible by the number of processors

//System dimensions are (maxi, maxj)
//Global array indices range from (0, 0) to (maxi-1, maxj-1)
//Internal points are from (1, 1) to (maxi-2, maxj-2)
//In each node mesh indices ranges from (node_mini, 0) to (node_maxi, maxj-1)

# include <math.h>
# include <stdlib.h>
# include <stdio.h>
#include <string.h>

#include <mpi.h>

void exchange_boundary_fields (void);
long double **matrix (long rows, long cols);
void initialize (void);
void printinfo (char *text);
void make_gnu_script(void);

//sides of the box
int imax, jmax;
//limits within which E and H are computed in each node
int startE, startH, endE, endHx, endHy;

//variables required for mpi related stuff
//size of spatial partition in mesh units
int block_size, node_maxi, node_mini;
int my_rank;
int nproc;

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

          //Write the Ez field at this step out to file
          //filename is created at runtime from n and my_rank
          sprintf (filename, "Ez%i.%i", n, my_rank);
          fp = fopen (filename, "w");
          for (i = node_mini; i <= node_maxi; i++) {
               for (j = 0; j < jmax; j++) {
                    fprintf (fp, "%i %i %Lf\n", i, j, ez[i][j]);
               }
               //blank line required by gnuplot
               fprintf(fp, " \n");
          }
          
          //close ez data file for this timestep             
          fclose (fp);
          fflush(fp);
          
          exchange_boundary_fields ();
          t = t + 1;
          
     }//for n - timestep ends
     
     if (my_rank == 0) {
          make_gnu_script ();
          system ("gnuplot -persist gnuem2-mpi.txt");
     }
     
     MPI_Finalize();
     return 0;
     
}

//end main

void exchange_boundary_fields (void)
{
int l;
MPI_Status status;
     
     //ez sent from process l+1 and received in process l
     for (l = 0; l < nproc - 1; l++) {
          if (my_rank == l) {
               MPI_Recv (&(ez[node_maxi+1][0]), jmax, MPI_LONG_DOUBLE, l+1, 0, MPI_COMM_WORLD, &status);
               
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
}

//end initialize

void make_gnu_script(void)
{
int i;
FILE *commandfile;
char text2[100];
     
     commandfile = fopen ("gnuem2-mpi.txt", "w");
     fprintf (commandfile, "set view 60, 30\n");
     
     for (n = 1; n <= nsteps; n++) {
          sprintf (text, "splot [][][-1:1]");
          for (i = 0; i < nproc; i++) {
               sprintf (text2, "'Ez%i.%i' w l, ", n, i);
               strcat (text, text2);
          }
          fprintf (commandfile, "%s\n", text);
     }
     
     fclose (commandfile);	
}

//end make_gnu_script

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