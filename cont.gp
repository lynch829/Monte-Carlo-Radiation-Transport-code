set terminal pngcairo size 1920,1080
set size square
set output 'contour.png'

#set view map
#unset surface
#set contours
#set cntrparam levels 10
#set cntrparam cubicspline

set xrange[46:54]
#set yrange[97:99]

set key outside

plot 'absorb.dat' using 2:1:3 matrix with image title ''
