set terminal pngcairo enhanced size 1920,1080
set size square
set output 'fluence.png'

set title 'Depth Resolved Fluence with Fresnel Reflections'
set ylabel 'Fluence[-]'
set xlabel 'z/cm'
set xtics ('0.0'99,'0.25'75,'0.50'50,'0.75'25,'1.0'0)
set logscale y
set ytics 1, 2, 10
set yrange[0.5:10]
set xrange [0:99] reverse

plot 'fres-same-sur.dat' w lp title 'n_1=curent', 'fres-diff-test.dat' w lp title 'n_1=1.37','fres-same-sur1.dat' w lp
