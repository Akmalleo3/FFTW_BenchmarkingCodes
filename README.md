# FFTW benchmarking codes

Code to time the roughly 100 us overhead in calling fftw plan when it
looks up stored wisdom to figure out it has already done that plan.

Example compile and demo (2000 repetitions, MEASURE, 3D, size 40^3, 4 threads):
```
export FINUFFT=locationofyourFINUFFTrepo
./compileFFTWBenchmark.sh
time ./fftwBenchmark 2000 1 3 40 4
```

Usage:
```
./fftwBenchMark [n_iterations] [FFTW_PLAN_TYPE] [DIMENSIONS] [SIZE (per dimension)]  [n_threads for fftw_execute (default omp_max_threads())] 

  plan types are {0:FFTW_MEASURE, 1:FFTW_ESTIMATE, 2: FFTW_PATIENT, 3: FFTW_EXHAUSTIVE}
```

### To do

* summarize results for various dims, plans, nthreads
