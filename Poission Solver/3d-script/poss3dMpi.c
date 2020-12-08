#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159275
#define MAXITS 1000
#define EPS 1.0e-6


void initialize (void);
void allocate (void);
void save_data (void);
double ***matrix3d (long rows, long cols, long depth);
void make_gnu_script(void);
void sor3d (double ***a, double ***b,
double ***c, double ***d,
double ***e, double ***f,
double ***g, double ***h,
double ***u,
int imax, int jmax, int kmax,
double rjac);
//coefficients
double ***a, ***b, ***c, ***d, ***e, ***g, ***h;
//potential
double ***u;
//charge
double ***f;
//side of box
int imax, jmax, kmax;
//spectral radius
double rjac;
int rank,processors;

int main (void)
{
  double t1=0;
  MPI_Init(0,0);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size( MPI_COMM_WORLD, &processors);

  initialize ();
  printf ("Iterating...\n");
  //solve Poisson's eqn
  t1 -= MPI_Wtime();
  sor3d (a, b, c, d, e, f, g, h, u, imax, jmax, kmax, rjac);
  t1 += MPI_Wtime();
  printf("Total time taken = %lf\n",t1);
  //what it says
  save_data ();
  printf ("Completed.\n");
  //make command file for gnuplot
  make_gnu_script();
  //show the results
  system ("gnuplot gnuscript3d.txt -persist");
  printf ("Displaying results.\n");
  return 0;
}

//end main
void initialize (void)
{
  int i, j, k;
  //system dimensyions
  imax = 100;
  jmax = 100;
  kmax = 100;
  //spectral radius of the Jacobi iteration matrix
  rjac = (cos (PI / imax) + cos (PI / jmax) + cos (PI / kmax))/3.0;
  //memory allocation clears values
  allocate ();
  //assign coefficients
  for (i = 0; i < imax; i++) {
    for (j = 0; j < jmax; j++){
      for (k = 0; k < kmax; k++){
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
  //create two parallel electrodes in the x-y-plane at k=0 and k=kmax
  for (i = 0; i < imax; i++) {
    for (j = 0; j < jmax; j++) {
      u[i][j][0] = -1.0e6;
      u[i][j][kmax] = 1.0e6;
    }
  }
}
//end intialize
void allocate (void)
{
  //coefficients
  a = matrix3d (imax, jmax, kmax);
  b = matrix3d (imax, jmax, kmax);
  c = matrix3d (imax, jmax, kmax);
  d = matrix3d (imax, jmax, kmax);
  f = matrix3d (imax, jmax, kmax);
  g = matrix3d (imax, jmax, kmax);
  //charge source
  h = matrix3d (imax, jmax, kmax);
  //central point
  e = matrix3d (imax, jmax, kmax);
  //potential
  u = matrix3d (imax, jmax, kmax);
}

//end allocate
void sor3d (double ***a, double ***b, double ***c,double ***d, double ***e, double ***f,double ***g, double ***h, double ***u,int imax, int jmax, int kmax,double rjac)
{
  int i, j, k, n, oddeven;
  double anorm, anormf=0.0, omega=1.0, resid;

  for (i = 0; i <= imax; i++)
    for (j = 0; j <= jmax; j++)
      for (k = 0; k <= kmax; k++)
        anormf += fabs(u[i][j][k])+fabs(h[i][j][k]);
  if (anormf == 0.0)
    return;
    for (n = 1; n <= MAXITS; n++) {
      //printf ("n=%i\n", n);
      oddeven = n % 2;
      anorm = 0.0;
      for (i = 1; i < imax; i++) {
        for (j = 1; j < jmax; j++) {
          for (k = 1; k < kmax; k++) {
            if (oddeven != ((i + j + k) % 2)) //skip this point
              continue;
            resid =a[i][j][k] * u[i+1][j][k] + b[i][j][k] * u[i-1][j][k] + c[i][j][k] * u[i][j+1][k] + d[i][j][k] * u[i][j-1][k] + f[i][j][k] * u[i][j][k+1] + g[i][j][k] * u[i][j][k-1] - e[i][j][k] * u[i][j][k] + h[i][j][k];
            anorm += fabs(resid);
            u[i][j][k] += omega * resid / e[i][j][k];
          }
        }
      }
      omega = 1.0/(1.0-0.25*rjac*rjac*omega);
//printf ("Iter %i: anorm = %f\n", n, anorm);
      if (anorm <= EPS * anormf) {
        printf ("Iterations = %i, Final anorm = %f, Initial anormf =%f\n", n, anorm, anormf);

      printf ("Final omega = %f\n", omega);
      return;
    }
  }
  printf ("MAXITS exceeded, anorm = %f anormf = %f\n", anorm, anormf);

}
//end sor3d

void save_data (void)
{
  FILE *output;
  int i, j;
  output = fopen ("solution3d-1.dat", "w");
  for (i = 0; i <= imax; i++) {
    for (j = 0; j <= jmax; j++) {
      fprintf (output, "%i %i %f\n", i, j, u[i][j][kmax/10]);
    }
    fprintf (output, "\n");
  }
  fclose (output);
  output = fopen ("solution3d-2.dat", "w");
  for (i = 0; i <= imax; i++) {
    for (j = 0; j <= jmax; j++) {
      fprintf (output, "%i %i %f\n", i, j, u[i][j][kmax/2]);
    }
    fprintf (output, "\n");
  }
  fclose (output);
  output = fopen ("solution3d-3.dat", "w");
  for (i = 0; i <= imax; i++) {
    for (j = 0; j <= jmax; j++) {
      fprintf (output, "%i %i %f\n", i, j, u[i][j][9*kmax/10]);
    }
    fprintf (output, "\n");
  }
  fclose (output);
}
//end save_data
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
//end matrix3d
void make_gnu_script(void)
{
  FILE *commandfile;
  commandfile = fopen ("gnuscript3d.txt", "w");
  fprintf (commandfile, "splot 'solution3d-1.dat' w l\n");
  fprintf (commandfile, "replot 'solution3d-2.dat' w l\n");

  fprintf (commandfile, "replot 'solution3d-3.dat' w l\n");
  fclose (commandfile);
}

// For Transfering Boundary Values
void boundaryTransfer(){
  for (int i = 0; i < processors-1; i++) {
    if(rank==i){

    }
    if(rank==i+1){
      
    }
  }
}
//end make_gnu_script
//end file
