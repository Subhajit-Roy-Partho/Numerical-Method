#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL,NULL);
	int rank, size, send_number=100,numbers[send_number];
	double mpi_time;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	for (int i = 0; i < send_number; ++i)
	{
		numbers[i] = rand()%1000;
	}
	srand(time(NULL));
	send_number = rand()%100;


	MPI_Barrier(MPI_COMM_WORLD);
	mpi_time -= MPI_Wtime();
	MPI_Bcast(numbers, send_number, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	mpi_time += MPI_Wtime();
	printf("Time = %lf for rank = %d\n",mpi_time,rank);
	MPI_Finalize();
	return 0;
}