#!/bin/bash

echo "------------Timing overhead of repeated construction of FFTW plans-----------------------"
echo ""
echo "A main objective of benchmarking the FFTW library was to determine the overhead"
echo "of repeated construction of fftw plans, in instances of repeated fft execution."
echo ""
echo "The fftwBenchmark compares such a scheme with the alternative: a single instantiation of a plan used"
echo "across repeated fft executions, using both the basic and advanced plan_many interface."
echo ""
echo ""


echo "-----------Plan overhead vs. dimension-------------------------------------"

D1=1
D2=2
D3=3
D1_size=900
D2_size=30
D3_size=10
n_iters=100

dimsize=(${D1_size} ${D2_size} ${D3_size})
for dimension in 1 2 3
do
    
    ./fftwBenchmark ${n_iters} 0 ${dimension} ${dimsize[${dimension}-1]} >> results.txt

    OUT="$(awk '/Mean time for later plan/ {print}' results.txt)"
    echo "${dimension} D: ${OUT}"
    rm results.txt
done


echo ""
echo ""
echo "-----------Plan overhead vs. thread count-------------------------------------"

for threadCount in 1 2 3 4 5 6 7 8
do
    ./fftwBenchmark ${n_iters} 0 ${D1} ${D1_size} ${threadCount} >> results.txt
    OUT="$(awk '/Mean time for later plan/ {print}' results.txt)"
    echo "${threadCount} threads, ${n_iters} trials, ${D1} D, size ${D1_size}: ${OUT}"
    rm results.txt
done

echo ""
echo ""
echo "-------Plan overhead vs. Number of FFTW Executions ----------------------------"

for trials in 10 50 100 500 1000 
do
    ./fftwBenchmark ${trials} 0 ${D1} ${D1_size}  >> results.txt
    OUT="$(awk '/Mean time for later plan/ {print}' results.txt)"
    echo "${trials} trials, ${D1} D, 8 threads, size ${D1_size}: ${OUT}"
    rm results.txt
done

echo ""
echo ""
echo "----Plan overhead vs. Plan Type {0: MEASURE, 1:ESTMATE, 2:PATIENT, 3:EXHAUSTIVE}--------"

for ptype in 0 1 #2 3
do
    ./fftwBenchmark ${n_iters} ${ptype} ${D2} ${D2_size}  >> results.txt
    OUT="$(awk '/Mean time for later plan/ {print}' results.txt)"
    echo "Type ${ptype}, ${n_iters} trials, ${D2} D, 8 threads, size ${D2_size}: ${OUT}"
    rm results.txt
done

echo "Types 2 and 3 were omitted due to unreasonable delays- comment back in at your own risk"



