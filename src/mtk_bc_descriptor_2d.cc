/*!
\file mtk_bc_descriptor_2d.cc

\brief Enforces boundary conditions in either the operator or the grid.

This class implements an interface for the user to specify boundary conditions
on 2D mimetic operators and the grids they are acting on.

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

#include "mtk_bc_descriptor_2d.h"

mtk::BCDescriptor2D::BCDescriptor2D():
  highest_order_diff_west(-1),
  highest_order_diff_east(-1),
  highest_order_diff_south(-1),
  highest_order_diff_north(-1),
  west_condition_(),
  east_condition_(),
  south_condition_(),
  north_condition_() {}

mtk::BCDescriptor2D::BCDescriptor2D(const mtk::BCDescriptor2D &desc) {}

mtk::BCDescriptor2D::~BCDescriptor2D() {}

void mtk::BCDescriptor2D::PushBackWestCoeff(mtk::CoefficientFunction2D cw) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cw == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_coefficients_.push_back(cw);

  highest_order_diff_west++;
}

void mtk::BCDescriptor2D::PushBackEastCoeff(mtk::CoefficientFunction2D ce) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(ce == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_coefficients_.push_back(ce);

  highest_order_diff_east++;
}

void mtk::BCDescriptor2D::PushBackSouthCoeff(mtk::CoefficientFunction2D cs) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cs == nullptr, __FILE__, __LINE__, __func__);
  #endif

  south_coefficients_.push_back(cs);

  highest_order_diff_south++;
}

void mtk::BCDescriptor2D::PushBackNorthCoeff(mtk::CoefficientFunction2D cn) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cn == nullptr, __FILE__, __LINE__, __func__);
  #endif

  north_coefficients_.push_back(cn);

  highest_order_diff_north++;
}

void mtk::BCDescriptor2D::set_west_condition_(
    mtk::Real (*west_condition)(mtk::Real xx, mtk::Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_condition_ = west_condition;
}

void mtk::BCDescriptor2D::set_east_condition_(
    mtk::Real (*east_condition)(mtk::Real xx, mtk::Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(east_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_condition_ = east_condition;
}

void mtk::BCDescriptor2D::set_south_condition_(
    mtk::Real (*south_condition)(mtk::Real xx, mtk::Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(south_condition == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  south_condition_ = south_condition;
}

void mtk::BCDescriptor2D::set_north_condition_(
    mtk::Real (*north_condition)(mtk::Real xx, mtk::Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(north_condition_ == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  north_condition_ = north_condition;
}

void mtk::BCDescriptor2D::ImposeOnLaplacianMatrix(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_cols() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Detect the order of accuracy of the Laplacian in the matrix.

  /// 2. Impose coefficients on the matrix.
}

void mtk::BCDescriptor2D::ImposeOnGrid(mtk::UniStgGrid2D &grid) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  #endif


}
