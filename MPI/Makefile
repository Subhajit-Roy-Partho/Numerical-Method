hello:
	mpicc hello.c -o hello;
	mpiexec --use-hwthread-cpus hello
clean:
	find . -type f -executable -exec rm '{}' \;