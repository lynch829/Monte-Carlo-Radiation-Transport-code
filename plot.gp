set term wxt size 1920,1080 persist
set size square
#load 'gnuplot-colorbrewer-master/diverging/Spectral.plt'
#set palette negative
set palette cubehelix
set xrange[0:200]
set yrange[0:200]
set key top right
datax=system("sed -n '13p' " . 'input.params' . "| awk '{print $1}'")
datay=system("sed -n '12p' " . 'input.params' . "| awk '{print $1}'")
datar=system("sed -n '14p' " . 'input.params' . "| awk '{print $1}'")
datap=system("sed -n '1p' " . 'numphotons.dat' . "| awk '{print $1}'")
set grid xtics ytics front lc rgb "white"
show grid
set title 'Detector at x='.datax.' y='.datay.' and a radius of '.datar.'cm detected '.datap.' photons'
numx=datax+0
numy=datay+0
rad = (datar*201)/2
#set cbrange[0:1]
#number after first is radius. formula for convo: (real_radius*nyg)/(2*xmax)
set object 1 circle at numx,numy size first rad fc rgb "white" front
#set pm3d map interpolate 0,0
plot 'density.dat' matrix with image
