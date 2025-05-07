#!/bin/bash

rm -f *.txt
make clean
make

for i in {1..30}
do 
    ./test_newton 100 >> newton_iter_100.txt
    sleep 1

    ./test_newton 1000 >> newton_iter_1000.txt
    sleep 1

    ./test_newton 10000 >> newton_iter_10000.txt
    sleep 1

    ./test_newton 100000 >> newton_iter_100000.txt
    sleep 1
done


gnuplot newton_time.gp

mkdir -p newton
mv *_100.txt *.png newton