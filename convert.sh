#!/bin/bash

rm -f *.txt
make clean
make

for i in {1..30}
do 
    ./test_convert >> fail.txt
    sleep 1
done

# plot
gnuplot conv_time.gp

# count failures
python3 count_f.py

# count errors
python3 count_err.py