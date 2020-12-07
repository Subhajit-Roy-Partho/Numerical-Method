#include <stdio.h>
#include <mpi.h>

void check_file (void);

int buf[10], rank, nproc;
char text[20];

int main(int argc, char *argv[])
{
	int i;
	MPI_File fh1, fh2;//file handles
	
	MPI_Status status;
	MPI_Offset offset;
	
	MPI_Init(0,0);//ignore (argc, argv)
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	
	for (i = 0; i < 10; i++){
		buf[i] = i * i + rank;
	}
	
	//create two files for writing
	sprintf (text, "test1.out");
	MPI_File_open(MPI_COMM_WORLD, text, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh1);
     sprintf (text, "test2.out");
     MPI_File_open(MPI_COMM_WORLD, text, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh2);
	
	if (rank == 1) {//only process 1 writes
		MPI_File_write(fh1, buf, 10, MPI_INT, MPI_STATUS_IGNORE);
	}
	
	offset = rank * 10 * sizeof (int);
	//all processes write to the same file
	MPI_File_write_at(fh2, offset, buf, 10, MPI_INT, &status);
	
	if (rank == 0)
		check_file ();//executed by process 0
	
	MPI_File_close(&fh1);
	MPI_File_close(&fh2);
	
	MPI_Finalize();
	
	return 0;
}

void check_file (void)
{
	int i, values[100];//this array should be large enough!
	FILE *fp;
	
	fp = fopen ("test1.out", "r");//this file written only by process 1
	fread (&values, sizeof (int), 10, fp);
	printf ("rank: %i: test1.out:\n", rank);
	for (i = 0; i < 10; i++){
		printf ("%i ", values[i]);
	}
	printf ("\n");
	fclose (fp);

	fp = fopen ("test2.out", "r");//this file was written by all processes
	fread (&values, sizeof (int), nproc*10, fp);
	printf ("rank: %i: test2.out:\n", rank);
	for (i = 0; i < nproc*10; i++){
		printf ("%i ", values[i]);
	}
	printf ("\n");
	fclose (fp);
	
}
