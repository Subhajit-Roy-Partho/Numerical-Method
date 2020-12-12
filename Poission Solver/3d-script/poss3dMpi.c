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
double ***matrix3d (long rows, long cols, long depth) //Create 3d matrix(x,y,z)
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

void allocate(void){//To allocate matrix using the above function
  a = matrix3d (dimx, dimy, dimz);//Cofficient for grid spacing adjustment
  b = matrix3d (dimx, dimy, dimz);
  c = matrix3d (dimx, dimy, dimz);
  d = matrix3d (dimx, dimy, dimz);
  f = matrix3d (dimx, dimy, dimz);
  g = matrix3d (dimx, dimy, dimz);
  //charge source
  h = matrix3d (dimx, dimy, dimz);
  //central point
  e = matrix3d (dimx, dimy, dimz);
  //potential
  u = matrix3d (dimx, dimy, dimz);
}

void initialise(void){
  for (int i = 1; i < dimx-1; i++) {
    for (int j = 0; j < dimy; j++){
      for (int k = 0; k < dimz; k++){
        //coefficients of the neighbouring points
        a[i][j][k] = 1.0;
        b[i][j][k] = 1.0;
        c[i][j][k] = 1.0;
        d[i][j][k] = 1.0;
        f[i][j][k] = 1.0;
        g[i][j][k] = 1.0;
        //coefficient of the central point
        e[i][j][k] = 6.0;
      }
    }
  }
  //create two parallel electrodes in the x-y-plane at k=0 and k=dimz
  for (int i = 1; i < dimx-1; i++) {
    for (int j = 0; j < dimy; j++) {
      u[i][j][0] = -1.0e6;
      u[i][j][dimz] = 1.0e6;
    }
  }
}



int startup(){
  if (dimx%processors == 0) {
    if(rank==0){//Printing only by a single processor
      printf("Dimensions are ok\n");
    }
    dimx= dimx/processors + 2; //Making the partation along the x axis with 2 extra values to hold the transfer data.
  }else{
    if(rank==1){//Printing only by a single processor
      printf("Dimensional Problem\n");
    }
    return 1;
  }
  allocate();
  initialise();
  return 0;
}


//Main Functions
int main() {

  MPI_Init(NULL,NULL); //Starting MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); // CPU rank is stored in rank
  MPI_Comm_size(MPI_COMM_WORLD,&processors); // Total number of CPUs.

  jrad = (cos (PI /dimx) + cos (PI /dimy) + cos (PI /dimz))/3.0;
  startup();


  MPI_Finalize(); //Stopping MPI
  return 0;
}
