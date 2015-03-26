set terminal pngcairo size 1920,1080
set size square
set output 'col.png'

set title '1 Pixel Slice Through Middle of Cube of Skin with Fresnel Reflections'
set ylabel 'Fluence Rate'
set xlabel 'z/cm'
set xtics ('0.0'200,'0.25'150,'0.50'100,'0.75'50,'1.0'0)
set xrange [] reverse

plot 'col1.dat' with l title 'no surface reflections','col2.dat' w lp title 'no fresnel','col3.dat' w l title 'surface reflections'
