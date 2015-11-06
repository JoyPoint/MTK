# gnuplot script for test suite: mtk_uni_stg_grid_1d.cc

#!/usr/bin/gnuplot

set terminal png
set output "mtk_uni_stg_grid_1d_test_04.png"
set title "Vector Field on a Uniform Staggered 1D grid"
set grid x
set key bmargin center horizontal
set xlabel "x"
set ylabel "v(x)"
plot 'mtk_uni_stg_grid_1d_test_04.dat' u 2:xticlabel(1) title "v(x)"
