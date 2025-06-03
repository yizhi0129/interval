#!/bin/bash

rm -f *.txt
make clean
make


for i in {1..10}
do 
    ./test_compression 53 >> mpfr_res_53.txt
    sleep 1
done

for i in {1..10}
do 
    ./test_compression 64 >> mpfr_res_64.txt
    sleep 1
done

for i in {1..10}
do 
    ./test_compression 128 >> mpfr_res_128.txt
    sleep 1
done

for i in {1..10}
do 
    ./test_compression 256 >> mpfr_res_256.txt
    sleep 1
done

for i in {1..10}
do 
    ./test_compression 512 >> mpfr_res_512.txt
    sleep 1
done

for i in {1..10}
do 
    ./test_compression2 >> compression2_res.txt
    sleep 1
done

gnuplot mpfr_time.gp

gnuplot multi_conv_time.gp

python3 mpfr_mem.py

python3 multi_mem.py

mv *.png multi_prec_convert/