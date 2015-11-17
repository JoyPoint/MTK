# \file mtk_uni_stg_grid_1d_test_04.gnu
#
# \brief gnuplot script for test suite mtk_uni_stg_grid_1d_test.cc
#
# Minimally-complete gnuplot script to visualize data files created by the
# mtk::UniStgGrid1D::WriteToFile method in the TestBindVectorField
# test implemented in the mtk_uni_stg_grid_1d_test.cc test suite.
#
# \warning Not intended to be a general solution gut a minimal guidance.
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

name = "mtk_uni_stg_grid_1d_test_04"

# wxt terminal (wxWidgets library) for live rendering.
set terminal wxt size 1024,768 enhanced font 'Verdana,10' persist

# png terminal for disk storage.
# set terminal png
# set output name.".png"

# Axes.
set xlabel "x"
set ylabel "u(x)"
set grid

# Title and legend.
set title "Vector Field on a Uniform Staggered 1D Grid"

set key bmargin center horizontal

plot 'mtk_uni_stg_grid_1d_test_04.dat' u 2:xticlabel(1) title "v(x)"
