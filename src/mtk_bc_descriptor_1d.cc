/*!
\file mtk_bc_descriptor_1d.cc

\brief Enforces boundary conditions in either the operator or the grid.

This class presents an interface for the user to specify boundary conditions
on 1D mimetic operators and the grids they are acting on.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
*/
/*
Copyright (C) 2015, Computational Science Research Center, San Diego State
University. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Modifications to source code should be reported to: esanchez@mail.sdsu.edu
and a copy of the modified files should be reported once modifications are
completed, unless these modifications are made through the project's GitHub
page: http://www.csrc.sdsu.edu/mtk. Documentation related to said modifications
should be developed and included in any deliverable.

2. Redistributions of source code must be done through direct
downloads from the project's GitHub page: http://www.csrc.sdsu.edu/mtk

3. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

4. Usage of the binary form on proprietary applications shall require explicit
prior written permission from the the copyright holders, and due credit should
be given to the copyright holders.

5. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

The copyright holders provide no reassurances that the source code provided does
not infringe any patent, copyright, or any other intellectual property rights of
third parties. The copyright holders disclaim any liability to any recipient for
claims brought against recipient by any third party for infringement of that
parties intellectual property rights.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "mtk_tools.h"

#include "mtk_bc_descriptor_1d.h"

void mtk::BCDescriptor1D::ImposeOnOperatorMatrix(mtk::DenseMatrix &matrix,
    const std::vector<mtk::Real> &west,
    const std::vector<mtk::Real> &east) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(west.size() > (unsigned int) matrix.num_cols(),
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east.size() > (unsigned int) matrix.num_cols(),
                      __FILE__, __LINE__, __func__);
  #endif

  /// 1. Assign the west array.

  for (unsigned int ii = 0; ii < west.size(); ++ii) {
    matrix.SetValue(0, ii, west[ii]);
  }

  /// 2. Assign the east array.

  for (unsigned int ii = 0; ii < east.size(); ++ii) {
    matrix.SetValue(matrix.num_rows() - 1,
                    matrix.num_cols() - 1 - ii,
                    east[ii]);
  }
}

void mtk::BCDescriptor1D::ImposeOnGrid(mtk::UniStgGrid1D &grid,
                                 const mtk::Real &omega,
                                 const mtk::Real &epsilon) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Assign the west condition.

  grid.discrete_field_u()[0] = omega;

  /// 2. Assign the east condition.

  grid.discrete_field_u()[grid.num_cells_x() + 2 - 1] = epsilon;
}
