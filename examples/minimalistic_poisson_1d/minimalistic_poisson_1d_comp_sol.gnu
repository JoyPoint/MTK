# gnuplot script for example: minimalistic_poisson_1d.cc

#!/usr/bin/gnuplot

set terminal png
set output "minimalistic_poisson_1d_comp_sol.png"
set title "Computed Solution"
set grid x
set key bmargin center horizontal
set xlabel "x"
set ylabel "s(x)"
plot 'minimalistic_poisson_1d_comp_sol.dat' u 2:xticlabel(1) title "s(x)"  w lp
