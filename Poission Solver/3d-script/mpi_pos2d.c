#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#define PI 3.14159275
//maximum no. of iterations
#define MAXITS 1000

void initialize (void);
void allocate (void);
void save_data (void);
double **matrix (long rows, long cols);
void make_gnu_script(void);
void sor2d (double **a, double **b,
		double **c, double **d,
		double **e, double **f,
		double **u,
		int imax, int jmax,
		double rjac);

void printinfo (char *text);
void exchange_boundary_phi (void);

//coefficients
double **a, **b, **c, **d, **e;
//potential
double **u;
//charge
double **f;

//sides of the box
int imax, jmax;

//Jacobi spectral radius
double rjac;

//storage for miscellaneous strings
char text[256];
FILE *output;

//variables required for mpi related stuff
//size of spatial partition in mesh units
int block_size;
int my_rank;
int nproc;

double t1, t2;

//these are the indices in each node w.r.t. the global mesh
int node_mini, node_maxi;

int main (int argc, char** argv)
{
     MPI_Init(&argc, &argv);
     MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
     MPI_Comm_size( MPI_COMM_WORLD, &nproc );

     //initialize the system dimensions
     imax = 200;
     jmax = 200;

	//allocate memory for the arrays
	allocate ();
	//initialize coefficients and electrode potentials
	initialize ();
	//solve the Poisson's equation
     t1 = MPI_Wtime();
	sor2d (a, b, c, d, e, f, u, imax, jmax, rjac);
     t2 = MPI_Wtime();
     sprintf (text, "Time in solver %f\n", t2-t1);
     printinfo (text);

     //write data to file
     save_data ();

	//create script file for gnuplot
	//comment this out if gnuplot is not available
     if (my_rank == 0) {
          make_gnu_script();
          //invoke system to plot the results
          system ("gnuplot gnuscript2d.txt -persist");
          t2 = MPI_Wtime();
          printf ("0: Cumulative wall clock time: %f\n", t2-t1);
     }

	MPI_Finalize();
	return 0;
}

//end main

void allocate (void)
{
     //coefficients
	a = matrix (imax, jmax);
	b = matrix (imax, jmax);
	c = matrix (imax, jmax);
	d = matrix (imax, jmax);
	e = matrix (imax, jmax);
     //charges
	f = matrix (imax, jmax);
	//potentials
	u = matrix (imax, jmax);
}

//end allocate

void initialize (void)
{
	int i, j;
	float p;

     //mesh is spatially partitioned along x-direction
     //maxi MUST BE divisible by nproc

     //this is the size of the partition in each node along x-direction
     //it is the no. of mesh points along x in each node
     block_size = imax / nproc;

     //use rank of this node to find the the minimum and maximum of the global index along x-direction
     //in each node mesh index ranges from (node_mini, 0) to (node_maxi, maxj-1)
     node_mini = block_size * my_rank;
     node_maxi = node_mini + block_size - 1;

     sprintf (text, "Block size %i Min i %i, Maxi %i", block_size, node_mini, node_maxi);
     printinfo (text);

     //spectral radius
     rjac = 0.5 * (cos (PI / imax) + cos(PI / jmax));

	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax; j++){
			//coefficients of the neighbouring points
			a[i][j] = 1.0;
			b[i][j] = 1.0;
			c[i][j] = 1.0;
			d[i][j] = 1.0;

			//coefficient of the central point
			e[i][j] = -4.0;

			//charge sources
			f[i][j] = 0.0;

			//potentials
			u[i][j] = 0.0;
		}
	}
	//f[100][100] = 1000.0;
     //set boundary potentials
	for (i = 20; i < imax-20; i++) {
          u[i][0] = 500.0;
          u[i][jmax-1] = -500.0;
     }
     //synchronize all nodes
     MPI_Barrier(MPI_COMM_WORLD);

}

//end intialize

void sor2d(double **a, double **b,
	    double **c, double **d,
	    double **e, double **f,
	    double **u,
	    int imax, int jmax,
	    double rjac)
{
	int i, i1, i2, j, n;
     int jp1, jm1;
	int oddeven;
	double anormmin, anormlocal, anormglobal=0.0;
     double omega=1.0, resid;

     for (i = node_mini; i <= node_maxi; i++)
          for (j = 0; j <= jmax-1; j++)
               anormlocal += fabs(u[i][j]) + fabs(f[i][j]);

     sprintf (text, "Starting local anorm = %f\n", anormlocal);
     printinfo (text);

     //find the global norm
     MPI_Allreduce (&anormlocal, &anormglobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

     //iterations terminate when global norm goes below minimum norm
     anormmin = anormglobal * 1.0e-3;
	//SOR factor
	omega = 1.0 / (1.0 - 0.5 * rjac * rjac);

     //general case for all nodes
     i1 = node_mini;
     i2 = node_maxi;

     //particular case for rank 0
     if (my_rank == 0) {
          i1 = 1;
     }
     //particular case for last node
     if (my_rank == nproc - 1) {
          i2 = node_maxi - 1;
     }
     //particular case for one process
     if (nproc == 1) {
          i2 = node_maxi - 1;
     }

     //iterations starts here
	for (n = 1; n <= MAXITS; n++) {
		anormlocal = 0.0;
		//even numbered iterations have this equal to zero
		oddeven = n % 2;
          //sweep begins
		for (j = 1; j < jmax-1; j++) {
               jp1 = j + 1;
               jm1 = j - 1;
			for (i = i1; i <= i2; i++) {
				//odd points are skipped during even iterations, and vice versa
				if (oddeven != ((i+j) % 2)) {
					continue;
				}
				resid=
				a[i][j]*u[i+1][j] +
				b[i][j]*u[i-1][j] +
				c[i][j]*u[i][jp1] +
				d[i][j]*u[i][jm1] +
				e[i][j]*u[i][j] +
				f[i][j];

                    anormlocal += fabs(resid);//this value is different in different nodes
				u[i][j] -= omega*resid/e[i][j];
			}
		}//half sweep ends

		exchange_boundary_phi ();
          //find the global norm
          MPI_Allreduce (&anormlocal, &anormglobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		//check for convergence
		if (anormglobal < anormmin) {
               sprintf (text, "Ending local anorm = %f\n", anormlocal);
               printinfo (text);
               return;//all nodes will return at the same iteration
		}
		//new value of omega for next half iteration
		omega = 1.0 / (1.0 - 0.25 * rjac * rjac * omega);
	}
}

//end sor

void exchange_boundary_phi (void)
{
     int j, l;
     MPI_Status status;

     //boundary layers exchanged with processor on right
     for (l = 0; l < nproc - 1; l++) {
          if (my_rank == l) {
               //l sends layer node_maxi to l+1
               MPI_Send (&(u[node_maxi][0]), jmax, MPI_DOUBLE, l+1, 0, MPI_COMM_WORLD);
               //l receives into layer node_maxi+1 from l+1
               MPI_Recv (&(u[node_maxi+1][0]), jmax, MPI_DOUBLE, l+1, 0, MPI_COMM_WORLD, &status);
          }

          if (my_rank == l+1) {
               //l+1 receives into layer node_mini-1 from l
               MPI_Recv (&(u[node_mini-1][0]), jmax, MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &status);
               //l+1 sends layer node_mini to l
               MPI_Send (&(u[node_mini][0]), jmax, MPI_DOUBLE, l, 0, MPI_COMM_WORLD);
          }
     }
}

//end exchange_boundary_phi

void save_data (void)
{
	int i, j;

     sprintf (text, "solution.%i", my_rank);
	output = fopen (text, "w");
	for (i = node_mini; i <= node_maxi; i++) {
		for (j = 0; j < jmax; j++) {
               fprintf (output, "%i %i %f\n", i, j, u[i][j]);
          }
          fprintf (output, "\n");
     }
     fclose (output);
}
//end save_data

double **matrix (long rows, long cols)
{
	int i;
	double **p;

	p = (double **) calloc (rows, sizeof(double *));
	for (i = 0; i < rows; i=i+1) {
		p[i] = calloc(cols, sizeof (double));
	}

	return p;

}

//end matrix

void make_gnu_script(void)
{
     int i;
	FILE *commandfile;

	commandfile = fopen ("gnuscript2d.txt", "w");
	fprintf (commandfile, "set view 70, 130\n");
	fprintf (commandfile, "splot 'solution.0' w l title \"Node 0\"\n");
     for (i = 1; i < nproc; i++) {
          fprintf (commandfile, "replot 'solution.%i' w l title \"Node %i\"\n", i, i);
     }

	fclose (commandfile);
}

//end make_gnu_script

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

//end file
