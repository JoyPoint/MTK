# \file 2d_poisson_known_comp_sol.gnu
#
# \brief gnuplot script for example 2d_poisson.cc
#
# Minimally-complete gnuplot script to visualize data files created by the
# mtk::UniStgGrid2D::WriteToFile method in the main module of the 2d_poisson.cc
# example.
#
# \warning Not intended to be a general solution but a minimal guidance.
#
# \sa https://github.com/ejspeiro/gnuplot-Scripts-Sci-Comp
#
# \author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

# Copyright (C) 2015, Computational Science Research Center, San Diego State
# University. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Modifications to source code should be reported to: esanchez@mail.sdsu.edu
# and a copy of the modified files should be reported once modifications are
# completed, unless these modifications are made through the project's GitHub
# page: http://www.csrc.sdsu.edu/mtk. Documentation related to said
# modifications should be developed and included in any deliverable.
#
# 2. Redistributions of source code must be done through direct
# downloads from the project's GitHub page: http://www.csrc.sdsu.edu/mtk
#
# 3. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 4. Usage of the binary form on proprietary applications shall require explicit
# prior written permission from the the copyright holders, and due credit should
# be given to the copyright holders.
#
# 5. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# The copyright holders provide no reassurances that the source code provided
# does not infringe any patent, copyright, or any other intellectual property
# rights of third parties. The copyright holders disclaim any liability to any
# recipient for claims brought against recipient by any third party for
# infringement of that parties intellectual property rights.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#!/usr/bin/gnuplot

reset

dat_file_name = "2d_poisson_known_comp_sol"
control_dat_file_name = "2d_poisson_known_sol"
computed_dat_file_name = "2d_poisson_comp_sol"

# Terminals.
# wxt terminal (wxWidgets library) for live rendering.
# set terminal wxt size 1024,768 enhanced font 'Verdana,10' persist

# png terminal for disk storage.
set terminal png
set output dat_file_name.".png"

# epslatex terminal for publication. (Proportions: 1024/768).
# set terminal epslatex standalone size 13cm,9.75cm color colortext 10
# set output dat_file_name.".tex"

set termoption dash

# Data visualization.
# View as a 2D map:
# set view map
# View as a 3D surface where z = u(x,y):
set view 60,340
# Style 1 for analytic/control data.
set style line 1 lt 2 lc rgb 'black' lw 1 pt 7 ps 0.5
# Style 2 for computed data.
set style line 2 lt 2 lc rgb 'black' lw 1 pt 7 ps 1
set palette defined (0 '#0000ff', 1 '#00ff00', 2 '#ff0000')

# Axes.
set autoscale fix
set grid
set format '$%g$'
set xlabel "$x$"
set x2tics
set ylabel "$y$"
set y2tics
set zlabel "$u(x,y)$"

# Title and legend.
set title "Control and Computed Solution"
set key bmargin center horizontal

# View coordinates of the centers:
splot control_dat_file_name.".dat" u 1:2:3:xticlabels(1):yticlabels(2) \
  w lp ls 1 palette, \
  computed_dat_file_name.".dat" u 1:2:3:xticlabels(1):yticlabels(2) \
  w p ls 2 palette

# Uncomment next line to view coordinates of the cell edges of the grid instead:
splot control_dat_file_name.".dat" u 1:2:3 w lp ls 1 palette, \
  computed_dat_file_name.".dat" u 1:2:3 w p ls 2 palette
