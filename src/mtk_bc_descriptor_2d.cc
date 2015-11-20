/*!
\file mtk_bc_descriptor_2d.cc

\brief Enforces boundary conditions in either the operator or the grid.

This class presents an interface for the user to specify boundary conditions
on 2D mimetic operators and the grids they are acting on.

<b>Def.</b> Let \f$ f \f$ be any scalar or vector field defined over a domain
\f$ \Omega \f$. We can specify any linear combination of \f$ f \f$ and its \f$
n \f$ derivatives to fulfill a condition, which we define as a **boundary
condition**:

\f[
\forall \mathbf{x} \in \partial\Omega:
  \sum_{i = 0}^{n}
    c_i(\mathbf{x})
        <\hat{\mathbf{n}}, \frac{\partial^i f}{\partial x^i}(\mathbf{x})> =
      \beta(\mathbf{x}).
\f]

This class receives information about the highest-order of differentiation,
\f$ n \f$, all possible coefficient functions, \f$ c_i(\mathbf{x}) \f$
for any subset of the boundary (south, north, west and east), and each condition
for any subset of the boundary, and takes care of assigning them to both, the
differentiation matrices and the grids.

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
  generate_space_(false),
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
  mtk::Tools::Prevent(highest_order_diff_west_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  west_coefficients_.push_back(cw);

  highest_order_diff_west_++;
}

void mtk::BCDescriptor2D::PushBackEastCoeff(mtk::CoefficientFunction2D ce) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(ce == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  east_coefficients_.push_back(ce);

  highest_order_diff_east_++;
}

void mtk::BCDescriptor2D::PushBackSouthCoeff(mtk::CoefficientFunction2D cs) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cs == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_south_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  south_coefficients_.push_back(cs);

  highest_order_diff_south_++;
}

void mtk::BCDescriptor2D::PushBackNorthCoeff(mtk::CoefficientFunction2D cn) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cn == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_north_ > 1,
                      __FILE__, __LINE__, __func__);
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

void mtk::BCDescriptor2D::ImposeOnSouthBoundary(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const int &order_accuracy) const {

  // At this point we have all of the information we need to fully impose the
  // south boundary condition:
  // 1. We have the collection of coefficients. The size of this collection
  // tells us the type of BC for this boundary.
  // 2. We have the grid that we can use to evaluate the coefficients at.
  // 3. We have the matrix where to place them.

  // For now, we are sure that we will NOT have more than 2 coefficients per
  // boundary. That is, we only support Robin type FOR NOW.

  if (generate_space_) {

    /// 1. Impose the Dirichlet condition first.

    // For the south-west corner:
    auto cc = (south_coefficients_[0])(grid.west_bndy(), grid.south_bndy());

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
      matrix.num_cols() << " columns." << std::endl;
    std::cout << "Setting at " << 0 << ' ' << 0 << std::endl;
    #endif

    matrix.SetValue(0, 0, cc);

    // Compute first centers per dimension.
    auto first_center_x = grid.west_bndy() + grid.delta_x()/mtk::kTwo;

    // For each entry on the diagonal (south boundary):
    for (int ii = 0; ii < grid.num_cells_x(); ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real xx = first_center_x + ii*grid.delta_x();
      // Evaluate and assign the Dirichlet coefficient.
      cc = (south_coefficients_[0])(xx, grid.south_bndy());

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Setting at " << ii + 1 << ' ' << ii + 1 << std::endl;
      #endif

      matrix.SetValue(ii + 1, ii + 1, cc);
    }

    // For the south-east corner:
    cc = (south_coefficients_[0])(grid.east_bndy(), grid.south_bndy());

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << grid.num_cells_x() + 1 << ' ' <<
      grid.num_cells_x() + 1 << std::endl;
    #endif

    matrix.SetValue(grid.num_cells_x() + 1, grid.num_cells_x() + 1, cc);

    /// 2. Impose the Neumann condition second.

    /// \todo Impose the Neumann conditions on every pole, for every scenario.
  } else {

    /// 1. Impose the Dirichlet condition first.

    // For each entry on the diagonal:
    for (int ii = 0; ii < grid.num_cells_x() + 2; ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real xx{(grid.discrete_domain_x())[ii]};
      // Evaluate and assign the Dirichlet coefficient.
      mtk::Real cc = (south_coefficients_[0])(xx,grid.south_bndy());
      matrix.SetValue(ii, ii, cc);
    }

    /// 2. Impose the Neumann condition second.
  }
}

void mtk::BCDescriptor2D::ImposeOnNorthBoundary(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const int &order_accuracy) const {

  // At this point we have all of the information we need to fully impose the
  // north boundary condition:
  // 1. We have the collection of coefficients. The size of this collection
  // tells us the type of BC for this boundary.
  // 2. We have the grid that we can use to evaluate the coefficients at.
  // 3. We have the matrix where to place them.

  // For now, we are sure that we will NOT have more than 2 coefficients per
  // boundary. That is, we only support Robin type FOR NOW.

  int north_offset{(grid.num_cells_y() + 1)*(grid.num_cells_x() + 2)};

  if (generate_space_) {

    /// 1. Impose the Dirichlet condition first.

    // For the north-west corner:
    mtk::Real cc =
      (north_coefficients_[0])(grid.west_bndy(), grid.north_bndy());

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
      matrix.num_cols() << " columns." << std::endl;
    std::cout << "Setting at " << north_offset << ' ' << north_offset <<
      std::endl;
    #endif

    matrix.SetValue(north_offset, north_offset, cc);

    // Compute first centers per dimension.
    auto first_center_x = grid.west_bndy() + grid.delta_x()/mtk::kTwo;

    // For each entry on the diagonal (north boundary):
    for (int ii = 0; ii < grid.num_cells_x(); ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real xx = first_center_x + ii*grid.delta_x();
      // Evaluate and assign the Dirichlet coefficient.
      cc = (north_coefficients_[0])(xx, grid.north_bndy());

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Setting at " << north_offset + ii + 1 << ' ' <<
        north_offset + ii + 1 << std::endl;
      #endif

      matrix.SetValue(north_offset + ii + 1, north_offset + ii + 1, cc);
    }

    // For the north-east corner:
    cc = (north_coefficients_[0])(grid.east_bndy(), grid.north_bndy());

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << north_offset + grid.num_cells_x() + 1 <<
      ' ' << north_offset + grid.num_cells_x() + 1 << std::endl;
    #endif

    matrix.SetValue(north_offset + grid.num_cells_x() + 1,
                    north_offset + grid.num_cells_x() + 1, cc);

    /// 2. Impose the Neumann condition second.
  } else {

    /// 1. Impose the Dirichlet condition first.

    // For each entry on the diagonal:
    for (int ii = 0; ii < grid.num_cells_x() + 2; ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real xx{(grid.discrete_domain_x())[ii]};
      // Evaluate and assign the Dirichlet coefficient.
      mtk::Real cc = (north_coefficients_[0])(xx, grid.north_bndy());
      matrix.SetValue(north_offset + ii, north_offset + ii, cc);
    }

    /// 2. Impose the Neumann condition second.
  }
}

void mtk::BCDescriptor2D::ImposeOnWestBoundary(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const int &order_accuracy) const {

  // At this point we have all of the information we need to fully impose the
  // west boundary condition:
  // 1. We have the collection of coefficients. The size of this collection
  // tells us the type of BC for this boundary.
  // 2. We have the grid that we can use to evaluate the coefficients at.
  // 3. We have the matrix where to place them.

  if (generate_space_) {

    /// 1. Impose the Dirichlet condition first.

    // For the south-west corner:
    auto cc = (west_coefficients_[0])(grid.west_bndy(), grid.south_bndy());

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
      matrix.num_cols() << " columns." << std::endl;
    std::cout << "Setting at " << 0 << ' ' << 0 << std::endl;
    #endif

    /// \note As it can be seen, we must adopt a convention about how to treat
    /// the corners. Based on a reasoning with Otilio, we will take the
    /// arithmetic mean.
    matrix.SetValue(0, 0, (matrix.GetValue(0, 0) + cc)/mtk::kTwo);

    int west_offset{grid.num_cells_x() + 1};

    auto first_center_y = grid.south_bndy() + grid.delta_y()/mtk::kTwo;

    // For each west entry on the diagonal (west boundary):
    for (int ii = 0; ii < grid.num_cells_y(); ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real yy = first_center_y + ii*grid.delta_y();
      // Evaluate and assign the Dirichlet coefficient.
      cc = (west_coefficients_[0])(grid.west_bndy(), yy);

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Setting at " << west_offset + ii + 1 << ' ' <<
        west_offset + ii + 1 << std::endl;
      #endif

      matrix.SetValue(west_offset + ii + 1, west_offset + ii + 1, cc);

      west_offset += grid.num_cells_x() + 1;
    }

    // For the north-west corner:
    cc = (west_coefficients_[0])(grid.west_bndy(), grid.north_bndy());

    west_offset += grid.num_cells_x() + 1;
    int aux{west_offset};
    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << aux << ' ' << aux << std::endl;
    #endif
    matrix.SetValue(aux, aux, (matrix.GetValue(aux, aux) + cc)/mtk::kTwo);

    /// 2. Impose the Neumann condition second.
  } else {

    /// 1. Impose the Dirichlet condition first.

    int west_offset{grid.num_cells_x() + 1};
    // For each west entry on the diagonal:
    for (int ii = 0; ii < grid.num_cells_y() + 2; ++ii) {
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real yy{(grid.discrete_domain_y())[ii]};
      // Evaluate and assign the Dirichlet coefficient.
      mtk::Real cc = (west_coefficients_[0])(grid.west_bndy(), yy);
      mtk::Real aux =
        (matrix.GetValue(west_offset + ii, west_offset + ii) + cc)/mtk::kTwo;
      matrix.SetValue(west_offset + ii, west_offset + ii, aux);
      west_offset += grid.num_cells_x() + 1;
    }

    /// 2. Impose the Neumann condition second.
  }
}

void mtk::BCDescriptor2D::ImposeOnEastBoundary(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const int &order_accuracy) const {

  // At this point we have all of the information we need to fully impose the
  // east boundary condition:
  // 1. We have the collection of coefficients. The size of this collection
  // tells us the type of BC for this boundary.
  // 2. We have the grid that we can use to evaluate the coefficients at.
  // 3. We have the matrix where to place them.

  if (generate_space_) {

    /// 1. Impose the Dirichlet condition first.

    // For the south-east corner:
    auto cc = (east_coefficients_[0])(grid.east_bndy(), grid.south_bndy());

    int east_offset{grid.num_cells_x() + 1};
    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
      matrix.num_cols() << " columns." << std::endl;
    std::cout << "Setting at " << east_offset << ' ' << east_offset <<
      std::endl;
    #endif

    matrix.SetValue(east_offset,
                    east_offset,
                    (matrix.GetValue(east_offset, east_offset) + cc)/mtk::kTwo);

    auto first_center_y = grid.south_bndy() + grid.delta_y()/mtk::kTwo;

    // For each east entry on the diagonal (east boundary):
    for (int ii = 0; ii < grid.num_cells_y(); ++ii) {

      east_offset += grid.num_cells_x() + 1;

      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real yy = first_center_y + ii*grid.delta_y();
      // Evaluate and assign the Dirichlet coefficient.
      cc = (east_coefficients_[0])(grid.east_bndy(), yy);

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Setting at " << east_offset + ii + 1 << ' ' <<
        east_offset + ii + 1 << std::endl;
      #endif

      matrix.SetValue(east_offset + ii + 1, east_offset + ii + 1, cc);
    }

    // For the north-east corner:
    cc = (east_coefficients_[0])(grid.east_bndy(), grid.north_bndy());

    east_offset += grid.num_cells_x() + 1;
    east_offset += grid.num_cells_x() + 1;
    int aux{east_offset};
    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << aux << ' ' << aux << std::endl;
    #endif
    matrix.SetValue(aux, aux, (matrix.GetValue(aux, aux) + cc)/mtk::kTwo);

    /// 2. Impose the Neumann condition second.

  } else {

    /// 1. Impose the Dirichlet condition first.

    int east_offset{grid.num_cells_x() + 1};
    // For each west entry on the diagonal:
    for (int ii = 0; ii < grid.num_cells_y() + 2; ++ii) {
      east_offset += grid.num_cells_x() + 1;
      // Evaluate next set spatial coordinates to evaluate the coefficient.
      mtk::Real yy{(grid.discrete_domain_y())[ii]};
      // Evaluate and assign the Dirichlet coefficient.
      mtk::Real cc = (east_coefficients_[0])(grid.east_bndy(), yy);
      mtk::Real aux =
        (matrix.GetValue(east_offset + ii, east_offset + ii) + cc)/mtk::kTwo;
      matrix.SetValue(east_offset + ii, east_offset + ii, aux);
    }

    /// 2. Impose the Neumann condition second.

  }
}

void mtk::BCDescriptor2D::ImposeOnLaplacianMatrix(
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const int &order_accuracy) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(highest_order_diff_south_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_north_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_west_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.nature() != mtk::SCALAR,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_cols() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. If we have not bound anything to the grid, then we have to generate
  /// our collection of spatial coordinates, as we evaluate the coefficients.

  generate_space_ = !grid.Bound();

  /// 2. Assign values to implement south boundary condition.

  ImposeOnSouthBoundary(grid, matrix, order_accuracy);

  /// 3. Assign values to implement north boundary condition.

  ImposeOnNorthBoundary(grid, matrix, order_accuracy);

  /// 4. Assign values to implement west boundary condition.

  ImposeOnWestBoundary(grid, matrix, order_accuracy);

  /// 5. Assign values to implement east boundary condition.

  ImposeOnEastBoundary(grid, matrix, order_accuracy);
}

void mtk::BCDescriptor2D::ImposeOnGrid(mtk::UniStgGrid2D &grid) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  #endif


}
