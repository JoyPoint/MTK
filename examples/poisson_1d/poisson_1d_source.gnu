#!/usr/bin/gnuplot

set terminal png
set output "poisson_1d_source.png"
set title "Source Term"
set grid x
set key bmargin center horizontal
set xlabel "x"
set ylabel "s(x)"
plot 'poisson_1d_source.dat' u 2:xticlabel(1) title "s(x)" w lp
