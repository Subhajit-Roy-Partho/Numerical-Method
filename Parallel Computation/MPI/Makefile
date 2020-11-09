hello:
	mpicc hello.c -o hello;
	mpiexec --use-hwthread-cpus hello

send_receive:
	mpicc send_receive.c -o send_receive
	mpiexec --use-hwthread-cpus send_receive

clean:
	find . -type f -executable -exec rm '{}' \;