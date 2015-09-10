set terminal png
set output "poisson_1d_known_sol.png"
set title "Know Solution"
set grid x
set key bmargin center horizontal
set xlabel "x"
set ylabel "u(x)"
plot 'poisson_1d_known_sol.dat' u 2:xticlabel(1) title "u(x)"
