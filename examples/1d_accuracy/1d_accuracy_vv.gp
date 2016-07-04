# \file 1d_accuracy_vv.gp
#
# \brief Plots 1D vector fields that are used for control in numerical methods.
#
# When a 1D vector field is taken as control, use this script to plot it.
#
# \warning Not intended to be a general solution but a minimal guidance.
#
# \author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

# The MIT License (MIT)
#
# Copyright (c) 2015 Eduardo Sanchez
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#!/usr/bin/gnuplot

reset

dat_file_name = "1d_accuracy_vv"

# Terminals.

# Set wxt terminal (wxWidgets library) for live rendering.
set terminal wxt size 1024,768 enhanced font 'Verdana,10' persist

# Set png terminal for disk storage.
# set terminal png
# set output dat_file_name.".png"

# Set epslatex terminal for publication. (Proportions: 1024/768).
# set terminal epslatex standalone size 13cm,9.75cm color colortext 10
# set output dat_file_name.".tex"

set termoption dash

# Data manipulation.

scalex = 0.075

pl(xx) = scalex*(xx)/sqrt(xx**2)

# Data visualization.

set palette defined (0 '#0000ff', 1 '#00ff00', 2 '#ff0000')

# Style for analytic/control data on cell edges.
set style line 1 lt 2 lc rgb 'black' lw 1 pt 7 ps 0.5

# Axes.

set autoscale fix
set grid
set format "$%g$"
set xlabel "$x$"
set ylabel "$||v(x)||$"

# Title and Legend.

set title "Control 1D Vector Field"
unset key

# Plot!

plot \
  dat_file_name.".dat" u 1:2:2:xtic(1) w lp ls 1 palette
#   dat_file_name.".dat" u 1:2:(pl($3)):(pl($3)):($3) \
#     w vectors head filled palette title "$v(x,y)$"
