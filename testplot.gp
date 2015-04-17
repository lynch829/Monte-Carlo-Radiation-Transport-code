set terminal pngcairo size 1920,1080
set palette cubehelix
set size square
set cbrange[0:.004]
unset key
do for [i=1:20] {
set output sprintf('ani%03.0f.png',i)
plot ''.i.'.dat' matrix with image
}
