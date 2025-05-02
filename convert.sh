#!/bin/bash

rm -f *.txt
make clean
make

for i in {1..30}
do 
    ./test_convert
    sleep 1
done

# plot
gnuplot conv_time.gp
gnuplot conv_dilat.gp

# count errors
python3 count_err.py