# gnuplot script for test suite: mtk_uni_stg_grid_1d.cc

#!/usr/bin/gnuplot

set terminal png
set output "mtk_uni_stg_grid_1d_test_03.png"

set title "Scalar Field on a Uniform Staggered 1D Grid"
set key bmargin center horizontal

set grid x
set xlabel "x"
set ylabel "u(x)"
plot 'mtk_uni_stg_grid_1d_test_03.dat' u 2:xticlabel(1) title "u(x)"
