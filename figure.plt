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
plot "figure_data.dat" u 1:2 t 'Energy' w l ls 11, \
	"figure_data.dat" u 1:3 t 'Heat Capacity' w l ls 12, \
        "figure_data.dat" u 1:7 t 'Order Parameter 1' w l ls 13, \
        "figure_data.dat" u 1:8 t 'Order Parameter 2' w l ls 14, \
        "figure_data.dat" u 1:9 t 'Chi 1' w l ls 15, \
        "figure_data.dat" u 1:10 t 'Chi 2' w l ls 16

