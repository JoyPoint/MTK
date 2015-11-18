/*!
\file mtk_bc_descriptor_2d.h

\brief Imposes boundary conditions in either the operator or the grid.

This class presents an interface for the user to specify boundary conditions
on 2D mimetic operators and the grids they are acting on.

<b>Def.</b> Let \f$ f \f$ be any scalar or vector field defined over a domain
\f$ \Omega \f$. We can specify any linear combination of \f$ f \f$ and its \f$
n \f$ derivatives to fulfill a condition, which we define as a **boundary
condition**:

\f[
\forall \mathbf{x} \in \partial\Omega:
  \sum_{i = 0}^{n}
    c_i(\mathbf{x})\frac{\partial^i f}{\partial x^i}(\mathbf{x}) =
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

#ifndef MTK_INCLUDE_BC_DESCRIPTOR_2D_H_
#define MTK_INCLUDE_BC_DESCRIPTOR_2D_H_

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_2d.h"

namespace mtk{

/*!
\typedef CoefficientFunction2D

\ingroup c07-mim_ops

\brief A function of a BC coefficient evaluated on a 2D domain.
*/
typedef Real (*CoefficientFunction2D)(Real, Real);

class BCDescriptor2D {
 public:
  /// \brief Default constructor.
  BCDescriptor2D();

  /*!
  \brief Copy constructor.

  \param [in] desc Given 2D descriptor.
  */
  BCDescriptor2D(const BCDescriptor2D &desc);

  /// \brief Destructor.
  ~BCDescriptor2D();

  /*!
  \brief Push back coefficient function at west of lowest order diff. available.

  \param [in] cw Function \f$ c_w(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackWestCoeff(CoefficientFunction2D cw);

  /*!
  \brief Push back coefficient function at east of lowest order diff. available.

  \param [in] ce Function \f$ c_e(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackEastCoeff(CoefficientFunction2D ce);

  /*!
  \brief Push back coefficient function south of lowest order diff. available.

  \param [in] cs Function \f$ c_s(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackSouthCoeff(CoefficientFunction2D cs);

  /*!
  \brief Push back coefficient function north of lowest order diff. available.

  \param [in] cn Function \f$ c_n(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void PushBackNorthCoeff(CoefficientFunction2D cn);

  /*!
  \brief Set boundary condition at west.

  \param [in] west_condition \f$ \beta_w(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_west_condition_(Real (*west_condition)(Real xx, Real yy));

  /*!
  \brief Set boundary condition at east.

  \param [in] east_condition \f$ \beta_e(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_east_condition_(Real (*east_condition)(Real xx, Real yy));

  /*!
  \brief Set boundary condition at south.

  \param [in] south_condition \f$ \beta_s(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_south_condition_(Real (*south_condition)(Real xx, Real yy));

  /*!
  \brief Set boundary condition at north.

  \param [in] north_condition \f$ \beta_n(x,y):\Omega\mapsto\mathbb{R} \f$.
  */
  void set_north_condition_(Real (*north_condition)(Real xx, Real yy));

  /*!
  \brief Imposes the condition on the operator represented as matrix.

  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input Laplacian operator.
  */
  void ImposeOnLaplacianMatrix(const UniStgGrid2D &grid,
                               DenseMatrix &matrix) const;

  /*!
  \brief Imposes the condition on the grid.

  \param[in,out] grid Grid upon which impose the desired boundary condition.
  */
  void ImposeOnGrid(UniStgGrid2D &grid) const;

private:
  int highest_order_diff_west; ///< Highest order of differentiation for west.
  int highest_order_diff_east; ///< Highest order of differentiation for east.
  int highest_order_diff_south; ///< Highest order of differentiation for south.
  int highest_order_diff_north; ///< Highest order of differentiation for north.

  std::vector<CoefficientFunction2D> west_coefficients_;  ///< Coeffs. west.
  std::vector<CoefficientFunction2D> east_coefficients_;  ///< Coeffs. east.
  std::vector<CoefficientFunction2D> south_coefficients_; ///< Coeffs. south.
  std::vector<CoefficientFunction2D> north_coefficients_; ///< Coeffs. south.

  Real (*west_condition_)(Real xx, Real yy);   ///< Condition for west.
  Real (*east_condition_)(Real xx, Real yy);   ///< Condition for east.
  Real (*south_condition_)(Real xx, Real yy);  ///< Condition for south.
  Real (*north_condition_)(Real xx, Real yy);  ///< Condition for north.
};
}
#endif  // End of: MTK_INCLUDE_BC_DESCRIPTOR_2D_H_
