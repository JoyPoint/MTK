# gnuplot script for test suite: mtk_uni_stg_grid_2d.cc

#!/usr/bin/gnuplot

reset

set terminal png
set output "mtk_uni_stg_grid_2d_test_04.png"

set title "Scalar Field u(x,y) on a Uniform Staggered 2D Grid"
set key bmargin center horizontal

set grid x
set grid y
set xlabel "x"
set ylabel "y"

# Temporal solution by
# http://www.kleerekoper.co.uk/2014/05/how-to-create-heatmap-in-gnuplot.html
set view map
set dgrid3d
splot "mtk_uni_stg_grid_2d_test_04.dat" using 1:2:3 with pm3d
