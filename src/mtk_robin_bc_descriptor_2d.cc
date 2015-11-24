/*!
\file mtk_robin_bc_descriptor_2d.cc

\brief Impose Robin boundary conditions on the operators and on the grids.

This class presents an interface for the user to specify Robin boundary
conditions on 2D mimetic operators and the grids they are acting on.

<b>Def.</b> Let \f$ u(\mathbf{x},t):\Omega\times [t_0, t_n]\mapsto\mathbb{R} \f$
be the solution to an ordinary or partial differential equation of interest. We
say that \f$ u \f$ satisfies a **Robin boundary condition on**
\f$ \partial\Omega \f$ if and only if there exists
\f$ \beta(\mathbf{x},t):\Omega\times [t_0, t_n]\mapsto\mathbb{R} \f$ so that:
\f[
\forall t \in [t_0,t_n]\; \forall \mathbf{x} \in \partial\Omega:
  \delta(\mathbf{x},t)u(\mathbf{x},t) +
    \eta(\mathbf{x},t)(\hat{\mathbf{n}}\cdot\nabla u) = \beta(\mathbf{x},t).
\f]

Intuitively, a **Robin boundary condition** is a constraint that must be
satisfied by any linear combination of any scalar field \f$ u \f$ and its first
normal derivative, in order for \f$ u \f$ to represent a unique solution to a
given ordinary or partial differential equation of interest.

Instances of this class receive information about the coefficient functions
and each condition for any subset of the boundary (west, east, south and north
in 2D). These instances then handle the complexity of placing the coefficients
in the differentiation matrices and the conditions in the grids.

\sa http://mathworld.wolfram.com/NormalVector.html

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

#include "mtk_robin_bc_descriptor_2d.h"

mtk::RobinBCDescriptor2D::RobinBCDescriptor2D():
  highest_order_diff_west_(-1),
  highest_order_diff_east_(-1),
  highest_order_diff_south_(-1),
  highest_order_diff_north_(-1),
  west_condition_(),
  east_condition_(),
  south_condition_(),
  north_condition_() {}

mtk::RobinBCDescriptor2D::RobinBCDescriptor2D(
    const mtk::RobinBCDescriptor2D &desc):
  highest_order_diff_west_(desc.highest_order_diff_west_),
  highest_order_diff_east_(desc.highest_order_diff_east_),
  highest_order_diff_south_(desc.highest_order_diff_south_),
  highest_order_diff_north_(desc.highest_order_diff_north_),
  west_condition_(desc.west_condition_),
  east_condition_(desc.east_condition_),
  south_condition_(desc.south_condition_),
  north_condition_(desc.north_condition_) {}

mtk::RobinBCDescriptor2D::~RobinBCDescriptor2D() noexcept {}

int mtk::RobinBCDescriptor2D::highest_order_diff_west() const noexcept {

  return highest_order_diff_west_;
}

int mtk::RobinBCDescriptor2D::highest_order_diff_east() const noexcept {

  return highest_order_diff_east_;
}

int mtk::RobinBCDescriptor2D::highest_order_diff_south() const noexcept {

  return highest_order_diff_south_;
}

int mtk::RobinBCDescriptor2D::highest_order_diff_north() const noexcept {

  return highest_order_diff_north_;
}

void mtk::RobinBCDescriptor2D::PushBackWestCoeff(
    mtk::CoefficientFunction1D cw) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cw == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_west_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  west_coefficients_.push_back(cw);

  highest_order_diff_west_++;
}

void mtk::RobinBCDescriptor2D::PushBackEastCoeff(
    mtk::CoefficientFunction1D ce) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(ce == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  east_coefficients_.push_back(ce);

  highest_order_diff_east_++;
}

void mtk::RobinBCDescriptor2D::PushBackSouthCoeff(
    mtk::CoefficientFunction1D cs) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cs == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_south_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  south_coefficients_.push_back(cs);

  highest_order_diff_south_++;
}

void mtk::RobinBCDescriptor2D::PushBackNorthCoeff(
    mtk::CoefficientFunction1D cn) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(cn == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_north_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  north_coefficients_.push_back(cn);

  highest_order_diff_north_++;
}

void mtk::RobinBCDescriptor2D::set_west_condition(
    mtk::Real (*west_condition)(const mtk::Real &yy,
                                const mtk::Real &tt)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_condition_ = west_condition;
}

void mtk::RobinBCDescriptor2D::set_east_condition(
    mtk::Real (*east_condition)(const mtk::Real &yy,
                                const mtk::Real &tt)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(east_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_condition_ = east_condition;
}

void mtk::RobinBCDescriptor2D::set_south_condition(
    mtk::Real (*south_condition)(const mtk::Real &xx,
                                 const mtk::Real &tt)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(south_condition == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  south_condition_ = south_condition;
}

void mtk::RobinBCDescriptor2D::set_north_condition(
    mtk::Real (*north_condition)(const mtk::Real &xx,
                                 const mtk::Real &tt)) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(north_condition == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  north_condition_ = north_condition;
}

bool mtk::RobinBCDescriptor2D::ImposeOnSouthBoundaryNoSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  // For the south-west corner:
  auto cc = (south_coefficients_[0])(grid.west_bndy(), time);

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
    cc = (south_coefficients_[0])(xx, time);

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << ii + 1 << ' ' << ii + 1 << std::endl;
    #endif

    matrix.SetValue(ii + 1, ii + 1, cc);
  }

  // For the south-east corner:
  cc = (south_coefficients_[0])(grid.east_bndy(), time);

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Setting at " << grid.num_cells_x() + 1 << ' ' <<
    grid.num_cells_x() + 1 << std::endl;
  #endif

  matrix.SetValue(grid.num_cells_x() + 1, grid.num_cells_x() + 1, cc);

  if (highest_order_diff_south_ > 0) {

    /// 2. Impose the Neumann condition.
  }
  /// \todo Impose the Neumann conditions on every pole, for every scenario.

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnSouthBoundaryWithSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  /// \todo Impose Harmonic mean on the corners for the case when the generated
  /// space is available, for all poles.

  // For each entry on the diagonal:
  for (int ii = 0; ii < grid.num_cells_x() + 2; ++ii) {
    // Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real xx{(grid.discrete_domain_x())[ii]};
    // Evaluate and assign the Dirichlet coefficient.
    mtk::Real cc = (south_coefficients_[0])(xx, time);
    matrix.SetValue(ii, ii, cc);
  }

  if (highest_order_diff_south_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnNorthBoundaryNoSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  int north_offset{(grid.num_cells_y() + 1)*(grid.num_cells_x() + 2)};

  /// 1. Impose the Dirichlet condition first.

  // For the north-west corner:
  mtk::Real cc =
    (north_coefficients_[0])(grid.west_bndy(), time);

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
    cc = (north_coefficients_[0])(xx, time);

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << north_offset + ii + 1 << ' ' <<
      north_offset + ii + 1 << std::endl;
    #endif

    matrix.SetValue(north_offset + ii + 1, north_offset + ii + 1, cc);
  }

  // For the north-east corner:
  cc = (north_coefficients_[0])(grid.east_bndy(), time);

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Setting at " << north_offset + grid.num_cells_x() + 1 <<
    ' ' << north_offset + grid.num_cells_x() + 1 << std::endl;
  #endif

  matrix.SetValue(north_offset + grid.num_cells_x() + 1,
                  north_offset + grid.num_cells_x() + 1, cc);

  if (highest_order_diff_north_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnNorthBoundaryWithSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose Dirichlet condition.

  int north_offset{(grid.num_cells_y() + 1)*(grid.num_cells_x() + 2)};

  /// For each entry on the diagonal:
  for (int ii = 0; ii < grid.num_cells_x() + 2; ++ii) {
    /// Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real xx{(grid.discrete_domain_x())[ii]};
    /// Evaluate and assign the Dirichlet coefficient.
    mtk::Real cc = (north_coefficients_[0])(xx, time);
    matrix.SetValue(north_offset + ii, north_offset + ii, cc);
  }

  if (highest_order_diff_north_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnWestBoundaryNoSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  // For the south-west corner:
  auto cc = (west_coefficients_[0])(grid.south_bndy(), time);

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
    matrix.num_cols() << " columns." << std::endl;
  std::cout << "Setting at " << 0 << ' ' << 0 << std::endl;
  #endif

  /// \note As it can be seen, we must adopt a convention about how to treat
  /// the corners. Based on a reasoning with Otilio, we will take the
  /// **harmonic mean**.

  mtk::Real harmonic_mean = mtk::kOne/matrix.GetValue(0, 0) + mtk::kOne/cc;
  harmonic_mean = mtk::kTwo/harmonic_mean;

  matrix.SetValue(0, 0, harmonic_mean);

  int west_offset{grid.num_cells_x() + 1};

  auto first_center_y = grid.south_bndy() + grid.delta_y()/mtk::kTwo;

  // For each west entry on the diagonal (west boundary):
  for (int ii = 0; ii < grid.num_cells_y(); ++ii) {
    // Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real yy = first_center_y + ii*grid.delta_y();
    // Evaluate and assign the Dirichlet coefficient.
    cc = (west_coefficients_[0])(yy, time);

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << west_offset + ii + 1 << ' ' <<
      west_offset + ii + 1 << std::endl;
    #endif

    matrix.SetValue(west_offset + ii + 1, west_offset + ii + 1, cc);

    west_offset += grid.num_cells_x() + 1;
  }

  // For the north-west corner:
  cc = (west_coefficients_[0])(grid.north_bndy(), time);

  west_offset += grid.num_cells_x() + 1;
  int aux{west_offset};
  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Setting at " << aux << ' ' << aux << std::endl;
  #endif

  harmonic_mean = mtk::kOne/matrix.GetValue(aux, aux) + mtk::kOne/cc;
  harmonic_mean = mtk::kTwo/harmonic_mean;

  matrix.SetValue(aux, aux, harmonic_mean);

  if (highest_order_diff_west_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnWestBoundaryWithSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  int west_offset{grid.num_cells_x() + 1};
  // For each west entry on the diagonal:
  for (int ii = 0; ii < grid.num_cells_y() + 2; ++ii) {
    // Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real yy{(grid.discrete_domain_y())[ii]};
    // Evaluate and assign the Dirichlet coefficient.
    mtk::Real cc = (west_coefficients_[0])(yy, time);
    matrix.SetValue(west_offset + ii, west_offset + ii, cc);
    west_offset += grid.num_cells_x() + 1;
  }

  if (highest_order_diff_west_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnEastBoundaryNoSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  // For the south-east corner:
  auto cc = (east_coefficients_[0])(grid.south_bndy(), time);

  int east_offset{grid.num_cells_x() + 1};
  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Matrix has " << matrix.num_rows() << " rows and " <<
    matrix.num_cols() << " columns." << std::endl;
  std::cout << "Setting at " << east_offset << ' ' << east_offset <<
    std::endl;
  #endif

  mtk::Real harmonic_mean =
    mtk::kOne/matrix.GetValue(east_offset,east_offset) + mtk::kOne/cc;
  harmonic_mean = mtk::kTwo/harmonic_mean;

  matrix.SetValue(east_offset, east_offset, harmonic_mean);

  auto first_center_y = grid.south_bndy() + grid.delta_y()/mtk::kTwo;

  // For each east entry on the diagonal (east boundary):
  for (int ii = 0; ii < grid.num_cells_y(); ++ii) {

    east_offset += grid.num_cells_x() + 1;

    // Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real yy = first_center_y + ii*grid.delta_y();
    // Evaluate and assign the Dirichlet coefficient.
    cc = (east_coefficients_[0])(yy, time);

    #if MTK_DEBUG_LEVEL > 0
    std::cout << "Setting at " << east_offset + ii + 1 << ' ' <<
      east_offset + ii + 1 << std::endl;
    #endif

    matrix.SetValue(east_offset + ii + 1, east_offset + ii + 1, cc);
  }

  // For the north-east corner:
  cc = (east_coefficients_[0])(grid.north_bndy(), time);

  east_offset += grid.num_cells_x() + 1;
  east_offset += grid.num_cells_x() + 1;
  int aux{east_offset};
  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Setting at " << aux << ' ' << aux << std::endl;
  #endif

  harmonic_mean =
    mtk::kOne/matrix.GetValue(aux, aux) + mtk::kOne/cc;
  harmonic_mean = mtk::kTwo/harmonic_mean;

  matrix.SetValue(aux, aux, harmonic_mean);

  if (highest_order_diff_east_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnEastBoundaryWithSpace(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  /// 1. Impose the Dirichlet condition first.

  int east_offset{grid.num_cells_x() + 1};
  // For each west entry on the diagonal:
  for (int ii = 0; ii < grid.num_cells_y() + 2; ++ii) {
    east_offset += grid.num_cells_x() + 1;
    // Evaluate next set spatial coordinates to evaluate the coefficient.
    mtk::Real yy{(grid.discrete_domain_y())[ii]};
    // Evaluate and assign the arithmetic mean of Dirichlet coefficients.
    mtk::Real cc = (east_coefficients_[0])(yy, time);
    matrix.SetValue(east_offset + ii, east_offset + ii, cc);
  }

  if (highest_order_diff_east_ > 0) {

    /// 2. Impose the Neumann condition.
  }

  return true;
}

bool mtk::RobinBCDescriptor2D::ImposeOnLaplacianMatrix(
    const mtk::Lap2D &lap,
    const mtk::UniStgGrid2D &grid,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

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

  /// If we have not bound anything to the grid, then we have to generate our
  /// collection of spatial coordinates, as we evaluate the coefficients.

  bool success{true};

  if (!grid.Bound()) {
    success = ImposeOnSouthBoundaryNoSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnNorthBoundaryNoSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnWestBoundaryNoSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnEastBoundaryNoSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
  } else {
    success = ImposeOnSouthBoundaryWithSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnNorthBoundaryWithSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnWestBoundaryWithSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
    success = ImposeOnEastBoundaryWithSpace(lap, grid, matrix, time);
    #if MTK_DEBUG_LEVEL > 0
    if (!success) {
      return false;
    }
    #endif
  }

  return success;
}

void mtk::RobinBCDescriptor2D::ImposeOnGrid(
    mtk::UniStgGrid2D &grid,
    const mtk::Real &time) const {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.num_cells_y() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(west_condition_ == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_condition_ == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(south_condition_ == nullptr,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_condition_ == nullptr,
                      __FILE__, __LINE__, __func__);
  #endif

  /// 1. Impose assuming an scalar grid.
  if (grid.nature() == mtk::SCALAR) {

    /// 1.1. Impose south condition.

    /// 1.1.1. Impose south-west corner.
    mtk::Real xx = grid.west_bndy();
    (grid.discrete_field())[0] = south_condition_(xx, time);

    /// 1.1.2. Impose south border.
    xx = xx + grid.delta_x()/mtk::kTwo;
    // For every point on the south boundary:
    for (int ii = 0; ii < grid.num_cells_x(); ++ii) {
      (grid.discrete_field())[ii + 1] =
        south_condition_(xx + ii*grid.delta_x(), time);
    }

    /// 1.1.3. Impose south-east corner.
    xx = grid.east_bndy();
    (grid.discrete_field())[grid.num_cells_x() + 1] =
      south_condition_(xx, time);

    /// 1.2. Impose north condition.

    /// 1.2.1. Impose north-west corner.
    xx = grid.west_bndy();
    int north_offset{(grid.num_cells_y() + 1)*(grid.num_cells_x() + 2)};
    (grid.discrete_field())[north_offset] = north_condition_(xx, time);

    /// 1.2.2. Impose north border.
    xx = xx + grid.delta_x()/mtk::kTwo;
    for (int ii = 0; ii < grid.num_cells_x(); ++ii) {
      (grid.discrete_field())[north_offset + ii + 1] =
        north_condition_(xx + ii*grid.delta_x(), time);
    }

    /// 1.2.3. Impose north-east corner.
    xx = grid.east_bndy();
    (grid.discrete_field())[north_offset + grid.num_cells_x() + 1] =
        north_condition_(xx, time);

    /// 1.3. Impose west condition.

    /// 1.3.1. Impose south-west corner.
    /// \note As per discussion with Otilio, we will take the **arithmetic**
    /// **mean** of the values of the BCs at the corners.
    mtk::Real yy = grid.south_bndy();
    (grid.discrete_field())[0] =
      ((grid.discrete_field())[0] + west_condition_(yy, time))/mtk::kTwo;

    /// 1.3.2. Impose west border.
    int west_offset{grid.num_cells_x() + 1 + 1};
    yy = yy + grid.delta_y()/mtk::kTwo;
    for (int ii = 0; ii < grid.num_cells_y(); ++ii) {
      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Adding on " << west_offset << "-th position." << std::endl;
      #endif
      (grid.discrete_field())[west_offset] =
        west_condition_(yy + ii*grid.delta_y(), time);
      west_offset += grid.num_cells_x() + 1 + 1;
    }

    /// 1.3.3. Impose north-west corner.
    yy = grid.north_bndy();
    north_offset = (grid.num_cells_y() + 1)*(grid.num_cells_x() + 2);
    (grid.discrete_field())[north_offset] =
      ((grid.discrete_field())[north_offset] + west_condition_(yy, time))/
        mtk::kTwo;

    /// 1.4. Impose east condition.

    /// 1.4.1. Impose south-east corner.
    yy = grid.south_bndy();
    int east_offset{grid.num_cells_x() + 1};
    (grid.discrete_field())[east_offset] =
      ((grid.discrete_field())[east_offset] + east_condition_(yy, time))/
        mtk::kTwo;

    /// 1.4.2. Impose east border.
    yy = yy + grid.delta_y()/mtk::kTwo;
    for (int ii = 0; ii < grid.num_cells_y(); ++ii) {
      east_offset += grid.num_cells_x() + 1 + 1;
      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Adding on " << east_offset << "-th position." << std::endl;
      #endif
      (grid.discrete_field())[east_offset] =
        east_condition_(yy + ii*grid.delta_y(), time);
    }

    /// 1.4.3. Impose north-east corner.
    yy = grid.north_bndy();
    (grid.discrete_field())[north_offset + grid.num_cells_x() + 1] =
      ((grid.discrete_field())[north_offset + grid.num_cells_x() + 1] +
      east_condition_(yy, time))/mtk::kTwo;

  } else {

    /// 2. Impose assuming a vector grid.

    /// \todo Implement imposition for vector-valued grids. Need research here!
  }
}
