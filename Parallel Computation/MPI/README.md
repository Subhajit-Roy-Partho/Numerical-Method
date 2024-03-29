# MPI
#### Execution
```
mpicc hello.c -o hello; mpiexec --use-hwthread-cpus hello
```
<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->
- `--use-hwthread-cpus` is used for utilization of the Hyperthreaded cores.

#### Makefile

- `make hello` to compile and run hello.c program which will execute a simple code block in all the shared cores and would also give the rank and total number of cores in hyperthreaded system.

- `make clean` to clean all the executable files find should be installed.

- `make send_receive` compile and run send_receive.c. This is a tutorial on Dynamic query and probing.


### Functions

- `MPI_Init(NULL,NULL)` Start MPI
- `MPI_Comm_size(MPI_COMM_WORLD, &core_size)` total number of cores.
- `MPI_Comm_rank(MPI_COMM_WORLD, &core_rank)` serial number of the core.
- `MPI_Send(void* data,int count,MPI_Datatype datatype,int destination,int tag,MPI_Comm communicator)` example `MPI_Send(numbers,send_number,MPI_INT,1,0,MPI_COMM_WORLD)`
- `MPI_Status status` - Mpi status datatype.
- `int MPI_Get_count(const MPI_Status * status, MPI_Datatype datatype, int *count)` example `MPI_Get_count(&status,MPI_INT,&count)` - gives the count of the number of data packets received. &count stores the number of data-packets received.

- `int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
       MPI_Comm comm, MPI_Status * status)` example `MPI_Recv(numbers,100,MPI_INT,0,0,MPI_COMM_WORLD,&status)` - for receiving point to point send data.

- `int MPI_Barrier(MPI_Comm comm)` - barrier point for sync. :rage: Caution: If all the process doesn't hit the barrier at least once the program would freeze for  every and won't proceed.

- `int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)`
	- root is the rank of the sending unit. Encountered by others will receive the data send by the root process.

- `double MPI_Wtime(void)` returns current time in double.

- `int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
       void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)` 

- `MPI_Finalize()` End MPI

### Flags

`--oversubscribe`