#!/bin/bash

# User to set environment variable FINUFFT to point to that repo's main dir.
# Otherwise this default is used:
if [ -z "$FINUFFT" ]; then
    FINUFFT=~/work/finufft
fi

g++-8  -o fftwBenchmark fftwBenchmark.cpp $FINUFFT/src/utils.cpp -I /usr/include -I $FINUFFT/src/ -I /usr/lib/gcc/x86_64-linux-gnu/8/include -L /usr/lib/x86_64-linux-gnu -L $FINUFFT/lib/ -lfftw3_omp -lfftw3 -lfinufft -fopenmp -lm
