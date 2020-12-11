#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


#define PI 3.14159275 //Approximate value of PI
#define MAXITS 1000 //Maximum number of itteration
#define ERROR 1.0e-6 //Minimum norm value for stopping the itteration


double jrad; // Jacobi Spectral Radius
int rank,processors; //rank for storing the rank of the CPU, processors stores total number of processors unit available.
int dimx=200, dimy=200,dimz=200; // box size in terms of (x,y,z)
double ***a,***b,***c,***d,***e,***g,***h; //Cofficient for adjusting grid spacing for our case =1.
double ***u;// Main matrix containing all the potentials.
double ***f; //Provision for addition of extra charges in the system.
// Functions


int startup(){
  if (dimx%processors == 0) {

  }else{
    cout<<"Dimentional Problem\n";
  }
}


double ***matrix3d (long rows, long cols, long depth)
{
  int i, j;
  double ***m;
  //allocate the rows
  m = (double ***) calloc (rows+1, sizeof(double **));
  for (i = 0; i <= rows; i=i+1) {
    //loop over the rows and allocate the columns for each row
    m[i] = (double **) calloc(cols+1, sizeof (double*));
    //loop over each column and allocate the depth for each row-column
    for (j = 0; j <= cols; j++)
      m[i][j] = (double *) calloc (depth+1, sizeof (double));
  }
  return m;
}


//Main Functions
int main() {

  MPI_Init(NULL,NULL); //Starting MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); // CPU rank is stored in rank
  MPI_Comm_size(MPI_COMM_WORLD,&processors); // Total number of CPUs.

  jrad = (cos (PI /dimx) + cos (PI /dimy) + cos (PI /dimz))/3.0;



  MPI_Finalize(); //Stopping MPI
  return 0;
}
