# include <math.h>
# include <stdlib.h>
# include <stdio.h>
//length of the system
//change this to alter the system dimension
#define CELLS 200
int main (void)
{
//electric and magnetic fields
float Ex[CELLS], Hy[CELLS];
//variables
int n, k;
//centre of system, total no. of steps
int kc, nsteps;
//variable time and pulse height
float t, pulse;
//peak time of pulse and its spread
float t0, spread;
//data file

FILE *fp;
//variables for the application of boundary conditions
float ex_low_m1, ex_high_m1;
//cell width in metres and timestep in sec
float dx, dt;
//propagation factor for Mur ABC
float propfactvac;
//Center of the problem space
//change this to alter the initial position of the pulse
//EM pulse originates here
kc = CELLS / 2;
//timestep at which the pulse peaks
t0 = 40.0;
//Width of the incident pulse
//change this to alter the pulsewidth
spread = 12;
//the cell size is 1 cm
dx = 0.01;
//this timestep takes two steps to move the wave one cell width
dt = dx / (2.0 * 3.0e8);
//Initialize to 1 to avoid premature termination of the while loop
nsteps = 1;
//factor used in applying the Mur ABCs
propfactvac = -1.0 * (3.0e8 * dt - dx) / (3.0e8 * dt + dx);
while (nsteps > 0) {
//Initialize all fields to zero
for (k = 0; k < CELLS; k++) {
  Ex[k] = 0.0;
  Hy[k] = 0.0;
}
t = 0;
//boundary condition variables set to zero
ex_low_m1 = ex_high_m1 = 0.0;
//Ask the user how many steps to run the main FDTD loop
//Beware: input only integers, any other value leads to interesting

// results!

printf( "No. of steps (0 to quit, integers only--> ");
scanf("%d", &nsteps);
if (nsteps <= 0)
exit(1);
n = 0;
//open data file for writing; this overwrites any existing file
fp = fopen("ExHy", "w");
//Main FDTD Loop
for (n = 1; n <= nsteps; n++) {
//create a gaussian source
pulse = exp(-0.5 * (pow((t0 - t) / spread, 2.0)));
//t is the time parameter which varies the magnitude of the

// Gaussian

t = t + 1;
// Calculate the Ex field
//index k goes to from k=1 to k=CELLS-2 because
//at k=0 and k=CELLS-1 the field is assigned using Mur

// Absorbing Boundary Conditions

for (k = 1; k < CELLS - 1; k++) {
  Ex[k] = Ex[k] + 0.5 * (Hy[k-1] - Hy[k]);
}
//pulse adds to existing field, "soft source"
Ex[kc] = Ex[kc] + pulse;
//implementation of Mur Absorbing Boundary Conditions
//lower side
Ex[0] = ex_low_m1 - propfactvac * (Ex[1] - Ex[0]);
//store the next-to-boundary value at this time step in m1
ex_low_m1 = Ex[1];
//higher side
Ex[CELLS - 1] = ex_high_m1 - propfactvac * (Ex[CELLS - 2] -

Ex[CELLS - 1]);

//store the next-to-boundary value at this time step in m1
ex_high_m1 = Ex[CELLS - 2];
//end of implementation of Mur Absorbing Boundary Conditions
//Calculate the Hy field
for (k = 0; k < CELLS - 1; k++) {
  Hy[k] = Hy[k] + (3.0e8*dt/dx) * (Ex[k] - Ex[k + 1]);
}
//print out the Ex and Hy fields in the data file
for (k = 0; k < CELLS; k++) {
fprintf(fp, "%i %f %f\n", k, Ex[k], Hy[k]);
}
//blank line to separate blocks for gnuplot
fprintf(fp, "\n");
}//for n
//End of the Main FDTD Loop
fclose (fp);
fflush (fp);
// open the script file for gnuplot
fp = fopen ("gnuinp.txt", "w");
for (k = 1; k < nsteps; k++) {
//write the command to plot Ex
fprintf (fp, "plot [0:199][-2:2] 'ExHy' every 1:1:0:%i:%i:%i u 1:2 w l\n", k, CELLS-1, k);

}
//close script file
fclose (fp);
fflush(fp);
//invoke the system to execute gnuplot
system ("gnuplot gnuinp.txt -persist");
}//while loop ends
return 0;
}

//end main
//end file
