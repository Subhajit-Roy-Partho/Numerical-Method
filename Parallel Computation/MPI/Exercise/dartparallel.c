#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define sqr(x)	((x)*(x))
#define DARTS 5000
struct timespec t1, t2;

double dboard(int darts)
{
	double x_coord,
	y_coord,       
	pi,            
	r;             
	int score,     
	srand(time());
	score = 0;

	
	for (n = 1; n <= darts; n++) {
		
		r = 1.0 * rand () / RAND_MAX;
		x_coord = (2.0 * r) - 1.0;
		r = 1.0 * rand () / RAND_MAX;
		y_coord = (2.0 * r) - 1.0;
		if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
			score++;
	}
	pi = 4.0 * (double)score/(double)darts;
	return(pi);
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL,NULL);
	int size, rank,rounds= 100;
	double time=0,avepi=0;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	time -= MPI_Wtime();
	for(int i=0; i < rounds/size; i++){
		pi = dboard(DART);
		avepi = (avepi+pi)/2.0; 
	}
	if (rank != 0)
	{
		
	}

	time += MPI_Wtime();
	return 0;
}