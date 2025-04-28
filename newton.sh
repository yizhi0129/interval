#!/bin/bash

rm -f *.txt
make clean
make

for i in {1..30}
do 
    ./test_newton >> newton_iter.txt
    sleep 1
done