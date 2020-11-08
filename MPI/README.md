#### Execution
```
mpicc hello.c -o hello; mpiexec --use-hwthread-cpus hello
```
<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->
- `--use-hwthread-cpus` is used for utilization of the Hyperthreaded cores.

#### Makefile

- `make hello` to compile and run hello.c program which will execute a simple code block in all the shared cores and would also give the rank and total number of cores in hyperthreaded system.

- `make clean` to clean all the executable files find should be installed.


### Functions

- `MPI_Comm_size(MPI_COMM_WORLD, &core_size)` total number of cores.
- `MPI_Comm_rank(MPI_COMM_WORLD, &core_rank)` serial number of the core.
- `MPI_Send(void* data,int count,MPI_Datatype datatype,int destination,int tag,MPI_Comm communicator)` example 