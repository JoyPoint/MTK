# gnuplot script for example: poisson_1d.cc

#!/usr/bin/gnuplot

reset

set terminal png
set output "poisson_1d_comp_sol.png"
set title "Computed Solution"
set grid x
set key bmargin center horizontal
set xlabel "x"
set ylabel "s(x)"
plot 'poisson_1d_comp_sol.dat' u 2:xticlabel(1) title "s(x)"  w lp
