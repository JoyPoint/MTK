/*!
\file mtk_robin_bc_descriptor_2d.h

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

#ifndef MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_2D_H_
#define MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_2D_H_

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"
#include "mtk_lap_2d.h"
#include "mtk_uni_stg_grid_2d.h"

namespace mtk{

/*!
\typedef CoefficientFunction1D

\ingroup c07-mim_ops

\brief A function of a BC coefficient evaluated on a 1D domain and time.
*/
typedef Real (*CoefficientFunction1D)(const Real &xx, const Real &tt);

/*!
\class RobinBCDescriptor2D

\ingroup c07-mim_ops

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
*/
class RobinBCDescriptor2D {
 public:
  /*!
	\brief Default constructor.
	*/
  RobinBCDescriptor2D();

  /*!
  \brief Copy constructor.

  \param [in] desc Given 2D descriptor.
  */
  RobinBCDescriptor2D(const RobinBCDescriptor2D &desc);

  /*!
  \brief Destructor.
  */
  ~RobinBCDescriptor2D() noexcept;

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
  \brief Getter for the highest order of differentiation in the south boundary.

  \return Integer highest order of differentiation in the south boundary.
  */
  int highest_order_diff_south() const noexcept;

  /*!
  \brief Getter for the highest order of differentiation in the north boundary.

  \return Integer highest order of differentiation in the north boundary.
  */
  int highest_order_diff_north() const noexcept;

  /*!
  \brief Push back coefficient function at west of lowest order diff. available.

  \param [in] cw Coeff.
    \f$ c_w(y,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void PushBackWestCoeff(CoefficientFunction1D cw);

  /*!
  \brief Push back coefficient function at east of lowest order diff. available.

  \param [in] cw Coeff.
    \f$ c_e(y,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void PushBackEastCoeff(CoefficientFunction1D ce);

  /*!
  \brief Push back coefficient function south of lowest order diff. available.

  \param [in] cw Coeff.
    \f$ c_s(x,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void PushBackSouthCoeff(CoefficientFunction1D cs);

  /*!
  \brief Push back coefficient function north of lowest order diff. available.

  \param [in] cw Coeff.
    \f$ c_n(x,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void PushBackNorthCoeff(CoefficientFunction1D cn);

  /*!
  \brief Set boundary condition at west.

  \param [in] west_condition
    \f$ \beta_w(y,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void set_west_condition(Real (*west_condition)(const Real &yy,
                                                 const Real &tt)) noexcept;

  /*!
  \brief Set boundary condition at east.

  \param [in] east_condition
    \f$ \beta_e(y,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void set_east_condition(Real (*east_condition)(const Real &yy,
                                                 const Real &tt)) noexcept;

  /*!
  \brief Set boundary condition at south.

  \param [in] south_condition
    \f$ \beta_s(x,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void set_south_condition(Real (*south_condition)(const Real &xx,
                                                   const Real &tt)) noexcept;

  /*!
  \brief Set boundary condition at north.

  \param [in] north_condition
    \f$ \beta_n(x,t):\partial\Omega\times[t_0,t_n]\mapsto\mathbb{R} \f$.
  */
  void set_north_condition(Real (*north_condition)(const Real &xx,
                                                   const Real &tt)) noexcept;

  /*!
  \brief Imposes the condition on the operator represented as matrix.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnLaplacianMatrix(const Lap2D &lap,
                               const UniStgGrid2D &grid,
                               DenseMatrix &matrix,
                               const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the grid.

  \param[in,out] grid Grid upon which impose the desired boundary condition.
  \param[in] time Current time snapshot. Default is kZero.
  */
  void ImposeOnGrid(UniStgGrid2D &grid, const Real &time = kZero) const;

private:
  /*!
  \brief Imposes the condition on the south boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnSouthBoundaryNoSpace(const Lap2D &lap,
                                    const UniStgGrid2D &grid,
                                    DenseMatrix &matrix,
                                    const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the north boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnNorthBoundaryNoSpace(const Lap2D &lap,
                                    const UniStgGrid2D &grid,
                                    DenseMatrix &matrix,
                                    const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the west boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnWestBoundaryNoSpace(const Lap2D &lap,
                                   const UniStgGrid2D &grid,
                                   DenseMatrix &matrix,
                                   const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the east boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnEastBoundaryNoSpace(const Lap2D &lap,
                                   const UniStgGrid2D &grid,
                                   DenseMatrix &matrix,
                                   const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the south boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnSouthBoundaryWithSpace(const Lap2D &lap,
                                      const UniStgGrid2D &grid,
                                      DenseMatrix &matrix,
                                      const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the north boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnNorthBoundaryWithSpace(const Lap2D &lap,
                                      const UniStgGrid2D &grid,
                                      DenseMatrix &matrix,
                                      const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the west boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnWestBoundaryWithSpace(const Lap2D &lap,
                                     const UniStgGrid2D &grid,
                                     DenseMatrix &matrix,
                                     const Real &time = kZero) const;
  /*!
  \brief Imposes the condition on the east boundary.

  \param[in] lap Laplacian operator on the matrix.
  \param[in] grid Grid upon which impose the desired boundary condition.
  \param[in,out] matrix Input matrix with the Laplacian operator.
  \param[in] time Current time snapshot. Default is kZero.
  */
  bool ImposeOnEastBoundaryWithSpace(const Lap2D &lap,
                                     const UniStgGrid2D &grid,
                                     DenseMatrix &matrix,
                                     const Real &time = kZero) const;

  int highest_order_diff_west_;   ///< Highest order of differentiation west.
  int highest_order_diff_east_;   ///< Highest order of differentiation east.
  int highest_order_diff_south_;  ///< Highest order differentiation for south.
  int highest_order_diff_north_;  ///< Highest order differentiation for north.

  std::vector<CoefficientFunction1D> west_coefficients_;  ///< Coeffs. west.
  std::vector<CoefficientFunction1D> east_coefficients_;  ///< Coeffs. east.
  std::vector<CoefficientFunction1D> south_coefficients_; ///< Coeffs. south.
  std::vector<CoefficientFunction1D> north_coefficients_; ///< Coeffs. north.

  Real (*west_condition_)(const Real &xx, const Real &tt);  ///< Condition west.
  Real (*east_condition_)(const Real &xx, const Real &tt);  ///< Condition east.
  Real (*south_condition_)(const Real &yy, const Real &tt); ///< Cond. south.
  Real (*north_condition_)(const Real &yy, const Real &tt); ///< Cond. north.
};
}
#endif  // End of: MTK_INCLUDE_ROBIN_BC_DESCRIPTOR_2D_H_
