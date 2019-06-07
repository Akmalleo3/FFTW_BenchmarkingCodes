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

------------Timing overhead of repeated construction of FFTW plans-----------------------

A main objective of benchmarking the FFTW library was to determine the overhead
of repeated construction of fftw plans, in instances of repeated fft execution.

The fftwBenchmark compares such a scheme with the alternative: a single instantiation of a plan used
across repeated fft executions, using both the basic and advanced plan_many interface.


-----------Plan overhead vs. dimension-------------------------------------
1 D: Mean time for later plan (us): 114.836
2 D: Mean time for later plan (us): 247.011
3 D: Mean time for later plan (us): 201.449


-----------Plan overhead vs. thread count-------------------------------------
1 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 12.8457
2 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 32.8656
3 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 48.3508
4 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 58.1091
5 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 76.3638
6 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 67.4532
7 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 98.2477
8 threads, 100 trials, 1 D, size 900: Mean time for later plan (us): 112.168


-------Plan overhead vs. Number of FFTW Executions ----------------------------
10 trials, 1 D, 8 threads, size 900: Mean time for later plan (us): 108.48
50 trials, 1 D, 8 threads, size 900: Mean time for later plan (us): 111.833
100 trials, 1 D, 8 threads, size 900: Mean time for later plan (us): 111.83
500 trials, 1 D, 8 threads, size 900: Mean time for later plan (us): 114.97
1000 trials, 1 D, 8 threads, size 900: Mean time for later plan (us): 119.522


----Plan overhead vs. Plan Type {0: MEASURE, 1:ESTMATE, 2:PATIENT, 3:EXHAUSTIVE}--------
Type 0, 100 trials, 2 D, 8 threads, size 30: Mean time for later plan (us): 258.2
Type 1, 100 trials, 2 D, 8 threads, size 30: Mean time for later plan (us): 347.636
Types 2 and 3 were omitted due to unreasonable delays- comment back in at your own risk
