#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL,NULL);
	int core_size;
	MPI_Comm_size(MPI_COMM_WORLD, &core_size);
	int core_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &core_rank);
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	printf("Processor Name =%s, rank = %d, Total Size %d",processor_name,core_rank, core_size);
	printf("\n");
	MPI_Finalize();
	return 0;
}
