# \file mimetic_threshold_div_grad.gp
#
# \brief gnuplot script for example mimetic_threshold_div_grad.cc.
#
# \warning Not intended to be a general solution but a minimal guidance.
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

dat_file_name = "mimetic_threshold_div_grad"

# Terminals.
# wxt terminal (wxWidgets library) for live rendering.
# set terminal wxt size 1024,768 enhanced font 'Verdana,10' persist

# png terminal for disk storage.
# set terminal png
# set output dat_file_name.".png"

# epslatex terminal for publication. (Proportions: 1024/768).
# Moved to plot section... see lines 98 and ?? of this file...

set termoption dash

# Data visualization.
# set palette defined (0 '#ffffff', 1 '#00ff00', 2 '#ff0000')

# Styles for computed data.
set style line 1 lt 2 lc rgb 'black' lw 1 pt 5 ps 0.9
set style line 2 lt 2 lc rgb 'black' lw 1 pt 3 ps 0.9
set style line 3 lt 2 lc rgb 'black' lw 1 pt 65 ps 0.9
set style line 4 lt 2 lc rgb 'black' lw 1 pt 11 ps 0.9

# Axes.

# Uncomment to plot on log scale for each axis respectively:
set logscale x
# set logscale y
# set logscale xy

# set autoscale fix
set grid
set format '$%g$'
set xlabel "$\\epsilon$"
set ylabel "Number of feasible solutions"

# Title and legend.
set title "Number of feasible solutions as a function of $\\epsilon$ for $\\breve{\\mathbf{D}}^k_x$"
set key below box height 1

set terminal epslatex standalone size 13cm,9.75cm color colortext 10
set output dat_file_name."-div.tex"

plot \
  'mimetic_threshold_div_8.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 1  \
  title "$\\breve{\\mathbf{D}}^8_x$", \
  'mimetic_threshold_div_10.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 2  \
  title "$\\breve{\\mathbf{D}}^{10}_x$", \
  'mimetic_threshold_div_12.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 3  \
  title "$\\breve{\\mathbf{D}}^{12}_x$", \
  'mimetic_threshold_div_14.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 4  \
  title "$\\breve{\\mathbf{D}}^{14}_x$"

set title "Number of feasible solutions as a function of $\\epsilon$ for $\\breve{\\mathbf{G}}^k_x$"
set output dat_file_name."-grad.tex"

plot \
  'mimetic_threshold_grad_10.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 2  \
  title "$\\breve{\\mathbf{D}}^{10}_x$", \
  'mimetic_threshold_grad_12.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 3  \
  title "$\\breve{\\mathbf{D}}^{12}_x$", \
  'mimetic_threshold_grad_14.dat' u 1:2:2:xtic(1):ytic(2) w lp ls 4  \
  title "$\\breve{\\mathbf{D}}^{14}_x$"
