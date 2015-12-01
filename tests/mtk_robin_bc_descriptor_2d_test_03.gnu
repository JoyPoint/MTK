# \file mtk_robin_bc_descriptor_2d_test_03.gnu
#
# \brief gnuplot script for test suite mtk_bc_descriptor_2d_test.cc
#
# Minimally-complete gnuplot script to visualize data files created by the
# mtk::UniStgGrid2D::WriteToFile method in the TestImposeOnGrid test
# implemented in the mtk_bc_descriptor_2d_test.cc test suite.
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

name = "mtk_robin_bc_descriptor_2d_test_03"

# wxt terminal (wxWidgets library) for live rendering.
set terminal wxt size 1024,768 enhanced font 'Verdana,10' persist

# png terminal for disk storage.
# set terminal png
# set output name.".png"

# Data visualization.
set grid
# Comment next line for a 3D surface plot where z = u(x,y):
# set view map

set palette defined (0 '#0000ff', \
                     1 '#00ff00', \
                     2 '#ff0000')

# Axes.
set xlabel "x"
set x2tics

set ylabel "y"

set autoscale fix

# Title and legend.
set title "Scalar Field on a Uniform Staggered 2D Grid"

set key bmargin center horizontal

# View coordinates of the centers:
splot name.".dat" u 1:2:3:xticlabels(1):yticlabels(2) w points pt 7 palette \
  title "u(x,y)"

# View coordinates of the cell edges:
splot name.".dat" u 1:2:3 w points pt 7 palette title "u(x,y)"
