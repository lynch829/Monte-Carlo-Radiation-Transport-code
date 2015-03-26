set terminal pngcairo size 1920,1080
set size square
set output 'acol.png'

set title '1 Pixel Slice Through Middle of Cube of Skin '
set ylabel 'Absorbance'
set xlabel 'z/cm'
set xtics ('0.0'200,'0.25'150,'0.50'100,'0.75'50,'1.0'0)
set xrange [] reverse

plot 'absorb.dat' with l title ''
#,'col2.dat' w l title 'no fresnel','col3.dat' w l title 'surface reflections'
