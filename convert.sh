#!/bin/bash

rm -f *.txt
make clean
make

for i in {1..30}
do 
    ./test_convert >> error.txt
done

# plot
gnuplot convert.gp