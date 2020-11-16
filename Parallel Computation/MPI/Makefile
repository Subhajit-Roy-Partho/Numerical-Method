hello:
	mpicc hello.c -o hello;
	mpiexec --use-hwthread-cpus hello

send_receive:
	mpicc send_receive.c -o send_receive
	mpiexec --use-hwthread-cpus send_receive

send_length:
	mpicc send_length.c -o send_length
	mpiexec --use-hwthread-cpus send_length

broadcast:
	mpicc broadcast.c -o broadcast;
	mpiexec --use-hwthread-cpus broadcast

scatter-gather:
	mpic++ scatter-gather.c -o scatter-gather;
	mpiexec --use-hwthread-cpus scatter-gather

# mpirun -n 5 --use-hwthread-cpus --oversubscribe row-mat

clean:
	find . -type f -executable -exec rm '{}' \;