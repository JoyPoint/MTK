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
  highest_order_diff_west_(-1),
  highest_order_diff_east_(-1),
  highest_order_diff_south_(-1),
  highest_order_diff_north_(-1),
  west_condition_(),
  east_condition_(),
  south_condition_(),
  north_condition_() {}

mtk::BCDescriptor2D::BCDescriptor2D(const mtk::BCDescriptor2D &desc) {}

mtk::BCDescriptor2D::~BCDescriptor2D() noexcept {}

int mtk::BCDescriptor2D::highest_order_diff_west() const noexcept {

  return highest_order_diff_west_;
}

int mtk::BCDescriptor2D::highest_order_diff_east() const noexcept {

  return highest_order_diff_east_;
}

int mtk::BCDescriptor2D::highest_order_diff_south() const noexcept {

  return highest_order_diff_south_;
}

int mtk::BCDescriptor2D::highest_order_diff_north() const noexcept {

  return highest_order_diff_north_;
}

void mtk::BCDescriptor2D::PushBackWestCoeff(mtk::CoefficientFunction2D cw) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cw == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_coefficients_.push_back(cw);

  highest_order_diff_west_++;
}

void mtk::BCDescriptor2D::PushBackEastCoeff(mtk::CoefficientFunction2D ce) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(ce == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_coefficients_.push_back(ce);

  highest_order_diff_east_++;
}

void mtk::BCDescriptor2D::PushBackSouthCoeff(mtk::CoefficientFunction2D cs) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cs == nullptr, __FILE__, __LINE__, __func__);
  #endif

  south_coefficients_.push_back(cs);

  highest_order_diff_south_++;
}

void mtk::BCDescriptor2D::PushBackNorthCoeff(mtk::CoefficientFunction2D cn) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cn == nullptr, __FILE__, __LINE__, __func__);
  #endif

  north_coefficients_.push_back(cn);

  highest_order_diff_north_++;
}

void mtk::BCDescriptor2D::set_west_condition(
    mtk::Real (*west_condition)(mtk::Real xx, mtk::Real yy)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_condition_ = west_condition;
}

void mtk::BCDescriptor2D::set_east_condition(
    mtk::Real (*east_condition)(mtk::Real xx, mtk::Real yy)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(east_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_condition_ = east_condition;
}

void mtk::BCDescriptor2D::set_south_condition(
    mtk::Real (*south_condition)(mtk::Real xx, mtk::Real yy)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(south_condition == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  south_condition_ = south_condition;
}

void mtk::BCDescriptor2D::set_north_condition(
    mtk::Real (*north_condition)(mtk::Real xx, mtk::Real yy)) noexcept {

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
  mtk::Tools::Prevent(grid.nature() != mtk::SCALAR,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_cols() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. If we have not bound anything to the grid, then we have to generate
  /// our collection of spatial coordinates, as we evaluate the coefficients.

  bool generate_space{false};

  if (!grid.Bound()) {
    generate_space = true;
  }

  /// 1. Impose the sum of the coefficients on the south boundary.

  if (generate_space) {
    // For the south-west corner:
    mtk::Real cc{};
    for (unsigned int jj = 0; jj < south_coefficients_.size(); ++jj) {
      // Evaluate the coefficient on the grid.
      cc += (south_coefficients_[jj])(grid.west_bndy(), grid.south_bndy());
    }
    matrix.SetValue(0, 0, cc);

    // Compute first centers per dimension.
    #ifdef MTK_PRECISION_DOUBLE
    auto first_center_x = grid.west_bndy() + grid.delta_x()/2.0;
    #else
    auto first_center_x = grid.west_bndy() + grid.delta_x()/2.0f;
    #endif

    #ifdef MTK_PRECISION_DOUBLE
    auto first_center_y = grid.south_bndy() + grid.delta_y()/2.0;
    #else
    auto first_center_y = grid.south_bndy() + grid.delta_y()/2.0f;
    #endif

    // For each entry on the diagonal (south boundary):
    for (int ii = 0; ii < grid.num_cells_x(); ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient at.
      mtk::Real xx = first_center_x + ii*grid.delta_x();
      mtk::Real yy = first_center_y + ii*grid.delta_y();

      cc = 0;
      // For each coefficient to sum:
      for (unsigned int jj = 0; jj < south_coefficients_.size(); ++jj) {
        // Evaluate the coefficient on the grid and accumulate.
        cc += (south_coefficients_[jj])(xx, yy);
      }
      matrix.SetValue(ii + 1, ii + 1, cc);
    }

    // For the south-east corner:
    cc = 0;
    for (unsigned int jj = 0; jj < south_coefficients_.size(); ++jj) {
      // Evaluate the coefficient on the grid.
      cc += (south_coefficients_[jj])(grid.east_bndy(), grid.south_bndy());
    }
    matrix.SetValue(grid.num_cells_x() + 1, grid.num_cells_x() + 1, cc);

  } else {
    // For each entry on the diagonal:
    for (int ii = 0; ii < grid.num_cells_x() + 2; ++ii) {
      // For each coefficient to sum:
      mtk::Real xx{(grid.discrete_domain_x())[ii]};
      mtk::Real yy{(grid.discrete_domain_y())[ii]};
      mtk::Real cc{};
      for (unsigned int jj = 0; jj < south_coefficients_.size(); ++jj) {
        // Evaluate the coefficient on the grid.
        cc += (south_coefficients_[jj])(xx,yy);
      }
      matrix.SetValue(ii, ii, cc);
    }
  }


  /// 2. Impose the sum of the coefficients on the north boundary.
}

void mtk::BCDescriptor2D::ImposeOnGrid(mtk::UniStgGrid2D &grid) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  #endif


}
