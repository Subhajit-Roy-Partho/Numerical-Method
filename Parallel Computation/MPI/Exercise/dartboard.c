//Calculate pi using Dartboard Monte Carlo algorithm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double dboard (int darts);

#define DARTS 5000   	/* number of throws at dartboard */

struct timespec t1, t2;

int main(int argc, char *argv[])
{
	double pi;          /* average of pi after "darts" is thrown */
	double avepi = 0;       /* average pi value for all iterations */
	int i, n, rounds;
	
	printf ("1st\n");
	rounds = atoi (argv[1]);
	printf ("2nd");
	printf("Starting pi calculation for %i rounds.\n", rounds);
	
	clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &t1);//start time
	
	for (i = 0; i < rounds; i++) {
		pi = dboard(DARTS);
		avepi = ((avepi * i) + pi)/(i + 1); 
	}
	
	clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &t2);//end time
	printf ("time = %lf\n", (t2.tv_sec+(t2.tv_nsec*1e-9) - t1.tv_sec+(t1.tv_nsec*1e-9)));//time taken
	
	printf("Average value of pi = %10.8lf. ", avepi);
	printf("Real value of PI: 3.1415926535897 \n");
	printf ("Difference is %lf\n", 3.1415926535897 - avepi);
}

/*****************************************************************************
 * dboard
 *****************************************************************************/
#define sqr(x)	((x)*(x))

double dboard(int darts)
{
	double x_coord,/* x coordinate, between -1 and 1  */
	y_coord,       /* y coordinate, between -1 and 1  */
	pi,            /* pi  */
	r;             /* random number scaled between 0 and 1  */
	int score,     /* number of darts that hit circle */
	n;
	
	score = 0;
	
	/***********************************************************************
	 * Throw darts at board.  Done by generating random numbers
	 * between 0 and 1 and converting them to values for x and y
	 * coordinates and then testing to see if they "land" in
	 * the circle."  If so, score is incremented.  After throwing the
	 * specified number of darts, pi is calculated.  The computed value
	 * of pi is returned as the value of this function, dboard.
	 ************************************************************************/
	
	for (n = 1; n <= darts; n++) {
		//generate random numbers for x and y coordinates
		r = 1.0 * rand () / RAND_MAX;//between 0 and 1
		x_coord = (2.0 * r) - 1.0;//shifted to interval -1 to 1
		r = 1.0 * rand () / RAND_MAX;//between 0 and 1
		y_coord = (2.0 * r) - 1.0;//shifted to interval -1 to 1
		
		//if dart lands in circle, increment score
		if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
			score++;
	}
	
	//calculate pi
	pi = 4.0 * (double)score/(double)darts;
	return(pi);
} 