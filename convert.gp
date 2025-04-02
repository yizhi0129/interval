set terminal png

set output "convert_time.png"

set key top right horizontal
set ylabel "time in ms"
set xlabel "methods"

set logscale y

set xrange [0:6]
set yrange [0.02:0.6]

set format x ""

stats 'convert_time.txt' using 1 name 'S1' nooutput
stats 'convert_time.txt' using 2 name 'S2' nooutput
stats 'convert_time.txt' using 3 name 'S3' nooutput
stats 'convert_time.txt' using 4 name 'S4' nooutput
stats 'convert_time.txt' using 5 name 'S5' nooutput
stats 'convert_time.txt' using 6 name 'S6' nooutput
stats 'convert_time.txt' using 7 name 'S7' nooutput
stats 'convert_time.txt' using 8 name 'S8' nooutput
stats 'convert_time.txt' using 9 name 'S9' nooutput
stats 'convert_time.txt' using 10 name 'S10' nooutput
stats 'convert_time.txt' using 11 name 'S11' nooutput


plot 'convert_time.txt' using (0):5 with points pt 6 ps 1.5 lc rgb 'yellow' title "read radius", \
     '' using (1):1 with points pt 6 ps 1.5 lc rgb '#42b983' title "CR 1", \
     '' using (2):2 with points pt 6 ps 1.5 lc rgb '#428bb9' title "CR 2", \
     '' using (3):3 with points pt 6 ps 1.5 lc rgb '#a8a9a8' title "CR 3", \
     '' using (4):4 with points pt 6 ps 1.5 lc rgb '#b468c7' title "CR 4", \
     '' using (5):10 with points pt 6 ps 1.5 lc rgb '#f5a945' title 'CR 5', \
     '' using (1.2):6 with points pt 6 ps 1.5 lc rgb '#359436' title "IS 1", \
     '' using (2.2):7 with points pt 6 ps 1.5 lc rgb '#4258b9' title "IS 2", \
     '' using (3.2):8 with points pt 6 ps 1.5 lc rgb '#414241' title "IS 3", \
     '' using (4.2):9 with points pt 6 ps 1.5 lc rgb '#813594' title "IS 4", \
     '' using (5.2):11 with points pt 6 ps 1.5 lc rgb '#794a0c' title 'IS 5', \
     '' using (1):(S1_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (2):(S2_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (3):(S3_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (4):(S4_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (0):(S5_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (1.2):(S6_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (2.2):(S7_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (3.2):(S8_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (4.2):(S9_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (5):(S10_median) with points pt 7 ps 1.5 lc rgb 'red' notitle, \
     '' using (5.2):(S11_median) with points pt 7 ps 1.5 lc rgb 'red' notitle