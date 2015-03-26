set terminal wxt persist

set logscale y
#set yrange[1:]
#set xrange[-30:30]

plot 'Tr.dat' w lp, 'rr.dat' w lp
