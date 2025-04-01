#!/bin/bash

make clean
make

for i in {1..20}
do 
    ./test_convert
done

# plot
gnuplot convert.gp