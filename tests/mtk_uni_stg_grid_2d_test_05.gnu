# gnuplot script for test suite: mtk_uni_stg_grid_2d.cc

#!/usr/bin/gnuplot

reset

# set terminal png
# set output "mtk_uni_stg_grid_2d_test_05.png"

set terminal wxt size 800,600 enhanced font 'Verdana,10' persist

set title "Vector Field v(x,y) on a Uniform Staggered 2D Grid"
set key bmargin center horizontal

set grid x
set grid y
set xlabel "x"
set ylabel "y"

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')

scalex = 0.075
scaley = 0.075

pl(xx,yy) = scalex*(xx + yy)/sqrt(xx**2 + yy**2)
ql(xx,yy) = scaley*(xx + yy)/sqrt(xx**2 + yy**2)

plot 'mtk_uni_stg_grid_2d_test_05.dat' w points pointsize 2,\
     'mtk_uni_stg_grid_2d_test_05.dat' u 1:2:(pl($3,$4)):(ql($3,$4)):($3 + $4) w vectors head filled palette
