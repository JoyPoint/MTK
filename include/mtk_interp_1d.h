/*!
\file mtk_interp_1d.h

\brief Includes the definition of the class Interp1D.

Definition of a class that implements a 1D interpolation operator.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\author: Johnny Corbino - jcorbino at mail dot sdsu dot edu
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

#ifndef MTK_INCLUDE_INTERP_1D_H_
#define MTK_INCLUDE_INTERP_1D_H_

#include <iostream>
#include <iomanip>

#include "glpk.h"

#include "mtk_roots.h"
#include "mtk_enums.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"

namespace mtk {

/*!
\class Interp1D

\ingroup c07-mim_ops

\brief Implements a 1D interpolation operator.

This class implements a 1D interpolation operator.
*/
class Interp1D {
 public:
  /// \brief Output stream operator for printing.
  friend std::ostream& operator <<(std::ostream& stream, Interp1D &in);

  /// \brief Default constructor.
  Interp1D();

  /*!
  \brief Copy constructor.

  \param [in] interp Given interpolation operator.
  */
  Interp1D(const Interp1D &interp);

  /// \brief Destructor.
  ~Interp1D();

  /*!
  \brief Factory method to build operator.

  \return Success of the solution.
  */
  bool ConstructInterp1D(int order_accuracy = kDefaultOrderAccuracy,
                         mtk::DirInterp dir = mtk::DirInterp::SCALAR_TO_VECTOR);

  /*!
  \brief Returns coefficients for the interior of the grid.

  \return Coefficients for the interior of the grid.
  */
  Real *coeffs_interior() const;

  /*!
  \brief Returns the operator as a dense matrix.

  \return The operator as a dense matrix.
  */
  DenseMatrix ReturnAsDenseMatrix(const UniStgGrid1D &grid) const;

 private:
  DirInterp dir_interp_;  ///< Direction of interpolation.

  int order_accuracy_;  ///< Order of numerical accuracy of the operator.

  Real *coeffs_interior_; ///< Interior stencil.
};
}
#endif  // End of: MTK_INCLUDE_INTERP_1D_H_
