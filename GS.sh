#!/bin/bash

rm -f *.txt
make clean
make


for i in {1..10}
do 
    for j in 10 100 1000 10000 100000
    do
        k=0 
        echo "Running test_GS with N = $j, precision = $k"
        ./test_GS "$j" "$k"
        sleep 1
    done
done

for i in {1..10}
do 
    for j in 10 100 1000 10000 100000
    do
        k=1 
        echo "Running test_GS with N = $j, precision = $k"
        ./test_GS "$j" "$k"
        sleep 1
    done
done



# plot

gnuplot gs_time.gp

gnuplot gs_err.gp