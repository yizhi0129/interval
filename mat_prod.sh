#!/bin/bash

rm -f *.txt
make clean
make 

mkdir matprod_adj

for i in {1..5}
do 
    ./test_functions 4000
    sleep 1

    ./test_functions 3000
    sleep 1

    ./test_functions 2000
    sleep 1

    ./test_functions 1000
    sleep 1

    ./test_functions 500
    sleep 1

    ./test_functions 100
    sleep 1

    ./test_functions 50
    sleep 1

    ./test_functions 10
    sleep 1
done

gnuplot matprod_time.gp

gnuplot matprod_dilat.gp

mv err1.txt err2.txt *.png matprod_adj

python3 dil_freq.py

mv freq1.txt freq2.txt freq_adj