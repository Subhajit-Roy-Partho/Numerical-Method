#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


#define PI 3.14159275
//maximum no. of iterations
#define MAXITS 1000 
//acceptable total residue
//iterations terminate when the total residue goes below this value
#define EPS 1.0e-3

void initialize (void);
void allocate (void);
void save_data (void);
double **matrix (long rows, long cols);
void make_gnu_script(void);
void sor (double **a, double **b,
		double **c, double **d,
		double **e, double **f,
		double **u,
		int imax, int jmax,
		double rjac);
int check_this_point (int i, int j);

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

int main (void)
{
	int rank, size;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//initialize the system dimensions

	imax = 200;
	jmax = 200;
	
	//spectral radius
	rjac = 0.5 * (cos (PI / imax) + cos(PI / jmax));
	
	//allocate memory for the arrays
	allocate ();
	//initialize coefficients and electrode potentials
	initialize ();
	//solve the Poisson's equation
	sor (a, b, c, d, e, f, u, imax, jmax, rjac);
	//write data to file
	save_data ();
	//create script file for gnuplot
	//comment this out if gnuplot is not available
	make_gnu_script();
	//invoke system to plot the results
	system ("gnuplot gnuscript2d.txt -persist");
	MPI_Finalize();
	return 0;
}

//end main

void allocate (void)
{
	a = matrix (imax, jmax);
	b = matrix (imax, jmax);
	c = matrix (imax, jmax);
	d = matrix (imax, jmax);
	e = matrix (imax, jmax);
	f = matrix (imax, jmax);
	
	u = matrix (imax, jmax);
}

//end allocate

void initialize (void)
{
	int i, j;
	float p;
	
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
	// f[100][100] = 1000.0;
	//--------------------
	//plane parallel electrodes held at constant potential
	//comment out this whole block to initialize a semi circular electrode
	//and make a similar change in function check_this_point below
	//changing the electrode dimensions must be followed by a similar change in function check_this_point
	// for (i = 20; i <= 180; i++) {
	// 	for (j = 50; j <= 55; j++) {
	// 		u[i][j] = 500.0;
	// 	}
	// 	for (j = 145; j <= 150; j++) {
	// 		u[i][j] = -500.0;
	// 	}
	// }
	//--------------------
	

	for (int i = 0; i < imax; ++i)
	{
		u[i][0] = 500.0;
		// u[i][jmax-1] = -500.0;
	}
	//--------------------
	/*
	 *     //semi-circular electrode centred at (imax/2,jmax/2)
	 *     //comment out this whole block to initialize a plane parallel electrode
	 *     //and make a similar change in function check_this_point below
	 *     //changing the electrode dimensions must be followed by a similar change in function check_this_point
	 *     for(i= imax/4; i<=imax/2; i++){
	 *          for(j= jmax/4; j<=3*jmax/4; j++)   {
	 *               p = sqrt ((i-imax/2.0) * (i-imax/2.0) + (j-jmax/2.0) * (j-jmax/2.0));
	 *               if(p >= imax/4.0 && p < imax/4.0 + 2.0)
	 *                    u[i][j]= 1000;
}
}
*/
	//--------------------
}

//end intialize

void sor(double **a, double **b,
	    double **c, double **d,
	    double **e, double **f,
	    double **u,
	    int imax, int jmax,
	    double rjac)
{
	int i, j, k, n;
	int oddeven;
	double anorm, anormf=0.0, omega=1.0, resid;
	
	for (i = 0; i <= imax; i++)
		for (j = 0; j <= jmax; j++)
			anormf += fabs(f[i][j] + u[i][j]);
		printf ("anormf at start = %f\n", anormf);
	//SOR factor
	omega = 1.0 / (1.0 - 0.5 * rjac * rjac);
	
	for (n = 1; n <= MAXITS; n++) {
		anorm = 0.0;
		//even numbered iterations have this equal to zero
		oddeven = n % 2;
		for (j = 1; j < jmax; j++) {
			for (i = 1; i < imax; i++) {
				//if the current point is an electrode point, skip it
				k = check_this_point (i, j);
				if (k == 1)
					continue;
				//odd points are skipped during even iterations, and vice versa
				if (oddeven != ((i+j) % 2)) {
					continue;
				}
				resid=a[i][j]*u[i+1][j] +
                         b[i][j]*u[i-1][j] +
                         c[i][j]*u[i][j+1] +
                         d[i][j]*u[i][j-1] +
                         e[i][j]*u[i][j] +
                         f[i][j];
				anorm += fabs(resid);
				u[i][j] -= omega*resid/e[i][j];
			}
		}
		omega = 1.0 / (1.0 - 0.25 * rjac * rjac * omega);
		
		//uncomment the line below to see the change in the SOR factor and the total residue
		//printf ("omega = %f anorm = %f\n", omega, anorm);
		
		//check for convergence
		if (anorm < EPS*anormf) {
			printf ("n = %f anorm = %f\n", n*0.5, anorm);
			return;
		}
	}
	//message in the event that the iterations fail to converge
	printf ("MAXITS exceeded: anorm = %f anormf = %f\n", anorm, anormf);
}

//end sor

int check_this_point (int i, int j)
{
	float p;
	
	//check if the point falls on or within plane parallel electrodes
	if (i >=20 && i <= 180) {
		if (j >= 50 && j <= 55)
			return 1;
		if (j >= 145 && j <= 150)
			return 1;
	}
	else
		return 0;
	
	/*
	 *     //check if point falls within circular electrode
	 *     p = sqrt ((i-imax/2.0) * (i-imax/2.0) + (j-jmax/2.0) * (j-jmax/2.0));
	 *     if(p >= imax/4.0 && p < imax/4.0 + 2.0)
	 *          return 1;
	 *     else
	 *          return 0;
	 */
}

void save_data (void)
{
	FILE *output;
	int i, j;
	
	output = fopen ("solution2d.dat", "w");
	for (i = 0; i <= imax; i++) {
		for (j = 0; j <= jmax; j++) {
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
	
	p = (double **) calloc (rows+1, sizeof(double *));
	for (i = 0; i <= rows; i=i+1) {
		p[i] = calloc(cols+1, sizeof (double));
	}
	
	return p;
	
}

//end matrix

void make_gnu_script(void)
{
	FILE *commandfile;
	
	commandfile = fopen ("gnuscript2d.txt", "w");
	fprintf (commandfile, "set hidden3d\n");
	fprintf (commandfile, "splot 'solution2d.dat' w l\n");
	fclose (commandfile);	
}

//end make_gnu_script

//end file
