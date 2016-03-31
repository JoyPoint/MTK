/*!
\file mtk_robin_bc_descriptor_1d.cc

\brief Impose Robin boundary conditions on the operators and on the grids.

This class presents an interface for the user to specify Robin boundary
conditions on 1D mimetic operators and the grids they are acting on.

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

In a 1D context (\f$ \partial\Omega = \{a,b\}\subset\mathbb{R} \f$), this
condition can be written as follows:
\f[
  \delta_a(a,t)u(a,t) - \eta_a(a,t)u^\prime(a,t) = \beta_a(a,t),
\f]
\f[
  \delta_b(b,t)u(b,t) + \eta_b(b,t)u^\prime(b,t) = \beta_b(b,t).
\f]

Instances of this class receive information about the coefficient functions
and each condition for any subset of the boundary (west and east, in 1D).
These instances then handle the complexity of placing the coefficients in the
differentiation matrices and the condition in the grids.

\sa http://mathworld.wolfram.com/NormalVector.html

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
*/
/*
Copyright (C) 2016, Computational Science Research Center, San Diego State
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
#include "mtk_grad_1d.h"
#include "mtk_robin_bc_descriptor_1d.h"

mtk::RobinBCDescriptor1D::RobinBCDescriptor1D():
  highest_order_diff_west_(-1),
  highest_order_diff_east_(-1),
  west_condition_(nullptr),
  east_condition_(nullptr) {}

mtk::RobinBCDescriptor1D::RobinBCDescriptor1D(
    const mtk::RobinBCDescriptor1D &desc):
  highest_order_diff_west_(desc.highest_order_diff_west_),
  highest_order_diff_east_(desc.highest_order_diff_east_),
  west_condition_(desc.west_condition_),
  east_condition_(desc.east_condition_) {}

mtk::RobinBCDescriptor1D::~RobinBCDescriptor1D() noexcept {}

int mtk::RobinBCDescriptor1D::highest_order_diff_west() const noexcept {

  return highest_order_diff_west_;
}

int mtk::RobinBCDescriptor1D::highest_order_diff_east() const noexcept {

  return highest_order_diff_east_;
}

void mtk::RobinBCDescriptor1D::PushBackWestCoeff(
    mtk::CoefficientFunction0D cw) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(cw == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_west_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  west_coefficients_.push_back(cw);

  highest_order_diff_west_++;
}

void mtk::RobinBCDescriptor1D::PushBackEastCoeff(
    mtk::CoefficientFunction0D ce) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(ce == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ > 1,
                      __FILE__, __LINE__, __func__);
  #endif

  east_coefficients_.push_back(ce);

  highest_order_diff_east_++;
}

void mtk::RobinBCDescriptor1D::set_west_condition(
    mtk::Real (*west_condition)(const mtk::Real &tt)) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(west_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  west_condition_ = west_condition;
}

void mtk::RobinBCDescriptor1D::set_east_condition(
    mtk::Real (*east_condition)(const mtk::Real &tt)) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(east_condition == nullptr, __FILE__, __LINE__, __func__);
  #endif

  east_condition_ = east_condition;
}

bool mtk::RobinBCDescriptor1D::ImposeOnDivergenceMatrix(
    const mtk::Div1D &div,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(highest_order_diff_west_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_cols() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Impose Dirichlet coefficients.

  /// 1.1. Impose Dirichlet condition at the west.
  matrix.SetValue(0, 0, (west_coefficients_[0])(time));

  /// 1.2. Impose Dirichlet condition at the east.
  matrix.SetValue(matrix.num_rows() - 1,
                  matrix.num_cols() - 1,
                  (east_coefficients_[0])(time));

  /// 2. Impose Neumann coefficients.

  if (highest_order_diff_west_ > 0) {}

  return true;
}

bool mtk::RobinBCDescriptor1D::ImposeOnLaplacianMatrix(
    const mtk::Lap1D &lap,
    mtk::DenseMatrix &matrix,
    const mtk::Real &time) const {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(highest_order_diff_west_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(highest_order_diff_east_ == -1,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(matrix.num_cols() == 0, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Impose Dirichlet coefficients.

  /// 1.1. Impose Dirichlet condition at the west.
  matrix.SetValue(0, 0, (west_coefficients_[0])(time));

  /// 1.2. Impose Dirichlet condition at the east.
  matrix.SetValue(matrix.num_rows() - 1,
                  matrix.num_cols() - 1,
                  (east_coefficients_[0])(time));

  /// 2. Impose Neumann coefficients.

  if (highest_order_diff_west_ > 0) {

    /// 2.1. Create a mimetic gradient to approximate the first derivative.
    mtk::Grad1D grad;
    if (!grad.ConstructGrad1D(lap.order_accuracy(),
                              lap.mimetic_threshold())) {
      return false;
    }

    /// 2.2. Extract the coefficients approximating the boundary.

    /// \warning Coefficients returned by the mim_bndy getter are dimensionless!
    /// Therefore we must scale them by delta_x (from the grid), before adding
    /// to the matrix! But this information is in the given lap!
    mtk::DenseMatrix coeffs(grad.mim_bndy());

    mtk::Real idx = mtk::kOne/lap.delta();

    /// 2.3. Impose Neumann condition at the west.
    for (int ii = 0; ii < coeffs.num_cols(); ++ii) {
      /// 2.3.1. Get gradient coefficient and scale it.
      mtk::Real aux{idx*coeffs.GetValue(0, ii)};
      /// 2.3.2. Multiply times the coefficient for this boundary, times the
      /// unit normal for this boundary.
      mtk::Real unit_normal{-mtk::kOne};
      aux *= unit_normal*(west_coefficients_[1])(time);
      /// 2.3.3. Set the final value summing it with what is on the matrix.
      matrix.SetValue(0, ii, matrix.GetValue(0, ii) + aux);
    }

    /// 2.4. Impose Neumann condition at the east.

    /// \warning The Coefficients returned by the mim_bndy getter are those
    /// intended for the west boundary. We must enforce the
    /// center-skew-symmetry of the resulting operator by permuting their
    /// location in the matrix, and changing their sign.

    for (int ii = 0; ii < coeffs.num_cols(); ++ii) {
      /// 2.4.1. Get gradient coefficient and scale it.
      mtk::Real aux{idx*coeffs.GetValue(0, ii)};
      /// 2.4.2. Multiply times the coefficient for this boundary, times the
      /// unit normal for this boundary, and change the sign to enforce
      /// center-skew-symmetry.
      mtk::Real unit_normal{mtk::kOne};
      aux *= -unit_normal*(east_coefficients_[1])(time);
      /// 2.4.3. Set the final value summing it with what is on the matrix.
      matrix.SetValue(matrix.num_rows() - 1,
                      matrix.num_rows() - 1 - ii,
                      matrix.GetValue(matrix.num_rows() - 1,
                                      matrix.num_rows() - 1 -ii) + aux);
    }
  }

  return true;
}

void mtk::RobinBCDescriptor1D::ImposeOnGrid(
    UniStgGrid1D &grid,
    const mtk::Real &time) const {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(west_condition_ == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_condition_ == nullptr, __FILE__, __LINE__, __func__);
  #endif

  (grid.discrete_field())[0] = west_condition_(time);
  (grid.discrete_field())[grid.num_cells_x() + 1] = east_condition_(time);
}
