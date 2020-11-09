#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL,NULL);
	int rank, size, send_number=20;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		MPI_Send(&send_number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}else if(rank == 1){
		MPI_Recv(&send_number, 1, MPI_INT, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("The received data = %d\n",send_number);
	}
	MPI_Finalize();
	return 0;
}