#!/bin/bash

rm -f *.txt
make clean
make


for i in {1..8}
do 
    for j in 100000 10000 1000 100 10
    do
        for k in 0 1 2
        do
            echo "Running test_GS_CSR with N = $j, precision = $k"
            ./test_GS_CSR "$j" "$k"
            sleep 1
        done
    done
done

for i in {1..8}
do 
    for j in 100000 10000 1000 100 10
    do
        for k in 53 64 128 256
        do
            echo "Running test_GS_CSR with N = $j, precision = $k"
            ./test_GS_CSR "$j" "$k"
            sleep 1
        done
    done
done


# plot

gnuplot gs_csr_time.gp

gnuplot GS_csr_mpfi_time.gp

gnuplot gs_csr_err.gp

gnuplot GS_csr_mpfi_err.gp

mv *.png GS_CSR
