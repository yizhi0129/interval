#!/bin/bash

rm -f *.txt
make clean
make


for i in {1..10}
do 
    for j in 10 100 1000 10000 100000
    do
        for k in 53 64 128 256
        do
            echo "Running test_GS with N = $j, precision = $k"
            ./test_GS "$j" "$k"
            sleep 1
        done
    done
done




# plot

#gnuplot gs_time.gp

#gnuplot gs_err.gp