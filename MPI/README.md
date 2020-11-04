#### Execution
```
mpicc hello.c -o hello; mpiexec --use-hwthread-cpus hello
```
<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->
- `--use-hwthread-cpus` is used for utilization of the Hyperthreaded cores.