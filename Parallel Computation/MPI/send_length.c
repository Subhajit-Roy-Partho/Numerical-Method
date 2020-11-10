#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL,NULL);
	int rank, size, send_number=100,numbers[send_number];
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		for (int i = 0; i < send_number; ++i)
		{
			numbers[i] = rand()%1000;
		}
		srand(time(NULL));
		send_number = rand()%100;
		MPI_Send(numbers,send_number,MPI_INT,1,0,MPI_COMM_WORLD);
		printf("Number of sent items = %d\n", send_number);
	}
	MPI_Finalize();
	return 0;
}