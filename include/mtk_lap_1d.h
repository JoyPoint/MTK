/*!
\file mtk_lap_1d.h

\brief Includes the definition of the class Lap1D.

This class implements a 1D Laplacian operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).

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

#ifndef MTK_INCLUDE_LAP_1D_H_
#define MTK_INCLUDE_LAP_1D_H_

#include <vector>

#include "mtk_dense_matrix.h"

#include "mtk_uni_stg_grid_1d.h"

namespace mtk {

/*!
\class Lap1D

\ingroup c07-mim_ops

\brief Implements a 1D mimetic Laplacian operator.

This class implements a 1D Laplacian operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).
*/
class Lap1D {
 public:
  /*!
	\brief Output stream operator for printing.
	*/
  friend std::ostream& operator <<(std::ostream& stream, Lap1D &in);

  /*!
  \brief Default constructor.
  */
  Lap1D();

  /*!
  \brief Copy constructor.

  \param [in] lap Given Laplacian.
  */
  Lap1D(const Lap1D &lap);

  /*!
  \brief Destructor.
  */
  ~Lap1D();

  /*!
  \brief Order of accuracy of the operator.

  \return Order of accuracy of the operator.
  */
  int order_accuracy() const;

  /*!
  \brief Mimetic threshold used in the CBS algorithm to construct this operator.

  \return Mimetic threshold used in the CBS algorithm to construct operator.
  */
  Real mimetic_threshold() const;

  /*!
  \brief Value of \f$ \Delta x \f$ used be scaled. If 0, then dimensionless.

  \return Value of \f$ \Delta x \f$ used be scaled. If 0, then dimensionless.
  */
  Real delta() const;

  /*!
  \brief Factory method implementing the CBS Algorithm to build operator.

  \return Success of the solution.
  */
  bool ConstructLap1D(int order_accuracy = kDefaultOrderAccuracy,
                      Real mimetic_threshold = kDefaultMimeticThreshold);

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
  \brief Return the operator as a dense array.

  \return The operator as a dense array.
  */
  const mtk::Real* data(const UniStgGrid1D &grid) const;

 private:
  int order_accuracy_;    ///< Order of numerical accuracy of the operator.
  int laplacian_length_;  ///< Length of the output array.

  Real *laplacian_;      ///< Output array containing the operator and weights.

  mutable Real delta_;     ///<< If 0.0, then this Laplacian is dimensionless.

  Real mimetic_threshold_;  ///<< Mimetic threshold.
  Real mimetic_measure_;    ///<< Mimetic measure.

  std::vector<Real> sums_rows_mim_bndy_; ///< Sum of each mimetic boundary row.
};
}
#endif  // End of: MTK_INCLUDE_LAP_1D_H_
