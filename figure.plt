#!/usr/bin/gnuplot
reset
  
set terminal pngcairo size 1000,400 enhanced font 'Verdana,10'
set output 'figure.png'

set border linewidth 1.5
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

set style fill transparent solid 0.5 noborder
set style function filledcurves y1=0
set clip two

load 'parula.pal'

set key outside right top
set lmargin 6
plot "test_output.dat" u 1:2 t 'Energy' w l ls 11, \
	"test_output.dat" u 1:3 t 'Heat Capacity' w l ls 12, \
	"test_output.dat" u 1:5 t 'Order Parameter' w l ls 13

