/*!
\file mtk_div_1d.h

\brief Declaration the definition of the class Div1D.

Definition of a class that implements a 1D divergence operator, constructed
using the Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).

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

#ifndef MTK_INCLUDE_DIV_1D_H_
#define MTK_INCLUDE_DIV_1D_H_

#include <iostream>
#include <iomanip>

#include <vector>

#include "glpk.h"

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"

namespace mtk {

/*!
\class Div1D

\ingroup c07-mim_ops

\brief Implements a 1D mimetic divergence operator.

This class implements a 1D divergence operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).
*/
class Div1D {
 public:
  /*!
	\brief Output stream operator for printing.
	*/
  friend std::ostream& operator <<(std::ostream& stream, Div1D &in);

  /*!
  \brief Default constructor.
  */
  Div1D();

  /*!
  \brief Copy constructor.

  \param [in] div Given divergence.
  */
  Div1D(const Div1D &div);

  /*!
  \brief Destructor.
  */
  ~Div1D();

  /*!
  \brief Factory method implementing the CBS Algorithm to build operator.

  \return Success of the construction.
  */
  bool ConstructDiv1D(int order_accuracy = kDefaultOrderAccuracy,
                      Real mimetic_threshold = kDefaultMimeticThreshold);

  /*!
  \brief Returns how many coefficients are approximating at the boundary.

  \return How many coefficients are approximating at the boundary.
  */
  int num_bndy_coeffs() const;

  /*!
  \brief Returns coefficients for the interior of the grid.

  \return Coefficients for the interior of the grid.
  */
  Real *coeffs_interior() const;

  /*!
  \brief Return collection of weights as computed by the CRSA.

  \return Collection of weights as computed by the CRSA.
  */
  Real *weights_crs(void) const;

  /*!
  \brief Return collection of weights as computed by the CBSA.

  \return Collection of weights as computed by the CBSA.
  */
  Real *weights_cbs(void) const;

  /*!
  \brief Return number of feasible solutions when using the CBSA for weights.

  \return Return number of feasible solutions when using the CBSA for weights.
  */
  int num_feasible_sols() const;

  /*!
  \brief Return collection of mimetic approximations at the boundary.

  \return Collection of mimetic approximations at the boundary.
  */
  DenseMatrix mim_bndy() const;

  /*!
  \brief Return collection of row-sums mimetic approximations at the boundary.

  \return Collection of row-sums mimetic approximations at the boundary.
  */
  std::vector<Real> sums_rows_mim_bndy() const;

  /*!
  \brief Returns mimetic measure of the operator.

  \return Real number which is the mimetic measure of the operator.
  */
  Real mimetic_measure() const;

  /*!
  \brief Return the operator as a dense matrix.

  \return The operator as a dense matrix.
  */
  DenseMatrix ReturnAsDenseMatrix(const UniStgGrid1D &grid) const;

  /*!
  \brief Returns the operator as a dimensionless dense matrix.

  \return The operator as a dimensionless dense matrix.
  */
  DenseMatrix ReturnAsDimensionlessDenseMatrix(int num_cells_x) const;

 private:
  /*!
  \brief Stage 1 of the CBS Algorithm.

  Compute the stencil approximating the interior of the staggered grid.
  */
  bool ComputeStencilInteriorGrid(void);

  /*!
  \brief Stage 2.1 of the CBS Algorithm.

  Compute a rational basis for the null-space of the Vandermonde matrix
  approximating at the west boundary.
  */
  bool ComputeRationalBasisNullSpace(void);

  /*!
  \brief Stage 2.2 of the CBS Algorithm.

  Compute the set of preliminary approximations on the boundary neighborhood.
  */
  bool ComputePreliminaryApproximations(void);

  /*!
  \brief Stage 2.3 of the CBS Algorithm.

  Compute the set of mimetic weights to impose the mimetic condition.
  */
  bool ComputeWeights(void);

  /*!
  \brief Stage 2.4 of the CBS Algorithm.

  Compute mimetic stencil approximating at boundary.
  */
  bool ComputeStencilBoundaryGrid(void);

  /*!
  \brief Stage 3 of the CBS Algorithm.

  Construct the output array with the operator and its weights.
  */
  bool AssembleOperator(void);

  int order_accuracy_;    ///< Order of numerical accuracy of the operator.
  int dim_null_;          ///< Dim. null-space for boundary approximations.
  int num_bndy_coeffs_;   ///< Req. coeffs. per bndy pt. uni. order accuracy.
  int divergence_length_; ///< Length of the output array.
  int minrow_;            ///< Row from the optimizer with the minimum rel. nor.
  int row_;               ///< Row currently processed by the optimizer.
  int num_feasible_sols_; ///< Number of feasible solutions for weights.

  DenseMatrix rat_basis_null_space_;  ///< Rational b. null-space w. bndy.

  Real *coeffs_interior_; ///< Interior stencil.
  Real *prem_apps_;       ///< 2D array of boundary preliminary approximations.
  Real *weights_crs_;     ///< Array containing weights from CRSA.
  Real *weights_cbs_;     ///< Array containing weights from CBSA.
  Real *mim_bndy_;        ///< Array containing mimetic boundary approximations.
  Real *divergence_;      ///< Output array containing the operator and weights.

  Real mimetic_threshold_;  ///<< Mimetic threshold.
  Real mimetic_measure_;    ///<< Mimetic measure.

  std::vector<Real> sums_rows_mim_bndy_; ///< Sum of each mimetic boundary row.
};
}
#endif  // End of: MTK_INCLUDE_DIV_1D_H_
