/*!
\file mtk_robin_bc_descriptor_1d.h

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

#include <vector>

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"
#include "mtk_div_1d.h"
#include "mtk_lap_1d.h"

#ifndef MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_1D_H_
#define MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_1D_H_

namespace mtk {
/*!
\typedef CoefficientFunction0D

\ingroup c07-mim_ops

\brief A function of a BC coefficient evaluated on a 0D domain and time.
*/
typedef Real (*CoefficientFunction0D)(const Real &tt,
																	    const std::vector<Real> &pp);

/*!
\class RobinBCDescriptor1D

\ingroup c07-mim_ops

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
differentiation matrices and the conditions in the grids.

\sa http://mathworld.wolfram.com/NormalVector.html
*/
class RobinBCDescriptor1D {
 public:
  /*!
	\brief Default constructor.
	*/
  RobinBCDescriptor1D();

  /*!
  \brief Copy constructor.

  \param [in] desc Given 1D descriptor.
  */
  RobinBCDescriptor1D(const RobinBCDescriptor1D &desc);

  /*!
  \brief Destructor.
  */
  ~RobinBCDescriptor1D() noexcept;

  /*!
  \brief Getter for the highest order of differentiation in the west boundary.

  \return Integer highest order of differentiation in the west boundary.
  */
  int highest_order_diff_west() const noexcept;

  /*!
  \brief Getter for the highest order of differentiation in the east boundary.

  \return Integer highest order of differentiation in the east boundary.
  */
  int highest_order_diff_east() const noexcept;

  /*!
  \brief Push back coefficient function at west of lowest order diff. available.

  \param [in] cw Function \f$ c_w(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackWestCoeff(CoefficientFunction0D cw);

  /*!
  \brief Push back coefficient function at east of lowest order diff. available.

  \param [in] ce Function \f$ c_e(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackEastCoeff(CoefficientFunction0D ce);

  /*!
  \brief Set boundary condition at west.

  \param [in] west_condition \f$ \beta_w(y,t):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_west_condition(Real (*west_condition)(const Real &tt)) noexcept;

  /*!
  \brief Set boundary condition at east.

  \param [in] east_condition \f$ \beta_e(y,t):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_east_condition(Real (*east_condition)(const Real &tt)) noexcept;

  /*!
  \brief Imposes the condition on the operator represented as matrix.

  \param[in] div Operator in the Matrix.
  \param[in,out] matrix Input Divergence operator.
  \param[in] time Current time snapshot. Default is kZero.

  \return Success of the imposition.
  */
  bool ImposeOnDivergenceMatrix(
  		const Div1D &div,
      DenseMatrix &matrix,
      const std::vector<Real> &parameters = std::vector<Real>(),
      const Real &time = mtk::kZero) const;

  /*!
  \brief Imposes the condition on the operator represented as matrix.

  \param[in] lap Operator in the Matrix.
  \param[in,out] matrix Input Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.

  \return Success of the imposition.
  */
  bool ImposeOnLaplacianMatrix(
  		const Lap1D &lap,
      DenseMatrix &matrix,
      const std::vector<Real> &parameters = std::vector<Real>(),
      const Real &time = mtk::kZero) const;
  /*!
  \brief Imposes the condition on the grid.

  \param[in,out] grid Grid upon which impose the desired boundary condition.
  \param[in] time Current time snapshot. Default is kZero.
  */
  void ImposeOnGrid(UniStgGrid1D &grid, const Real &time = mtk::kZero) const;

 private:
  int highest_order_diff_west_; 		///< Highest order of differentiation west.
  int highest_order_diff_east_; 		///< Highest order of differentiation east.

  std::vector<CoefficientFunction0D> west_coefficients_;  ///< Coeffs. west.
  std::vector<CoefficientFunction0D> east_coefficients_;  ///< Coeffs. east.

  Real (*west_condition_)(const Real &tt); ///< Condition for west.
  Real (*east_condition_)(const Real &tt); ///< Condition for east.
};
}
#endif  // End of: MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_1D_H_
