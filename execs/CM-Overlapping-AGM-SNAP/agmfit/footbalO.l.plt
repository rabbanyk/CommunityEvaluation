#
# N:115, E:613  (Mon Feb  9 18:19:26 2015)
#

set title "N:115, E:613 "
set key bottom right
set logscale x 10
set format x "10^{%L}"
set mxtics 10
set grid
set tics scale 2
set terminal png small size 1000,800
set output 'footbalO.l.png'
plot 	"footbalO.l.tab" using 5:6 title "C" with linespoints ,\
	"footbalO.l.tab" using 1:2 title "likelihood" with linespoints ,\
	"footbalO.l.tab" using 1:4 title "BIC" with linespoints ,\
	"footbalO.l.tab" using 1:3 title "Sigmoid (scaled)" with linespoints 
