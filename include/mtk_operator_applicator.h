/*!
\file mtk_operator_applicator.h

\brief Controls the process of applying a mimetic operator to a grid.



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

#ifndef MTK_INCLUDE_OPERATOR_APPLICATOR_H_
#define MTK_INCLUDE_OPERATOR_APPLICATOR_H_

#include <iostream>

#include "mtk_roots.h"
#include "mtk_enums.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"
#include "mtk_uni_stg_grid_2d.h"
#include "mtk_uni_stg_grid_3d.h"

namespace mtk {

/*!
\class OperatorApplicator

\ingroup c07-mim_ops

\brief Controls the process of applying a mimetic operator to a grid.


*/
class OperatorApplicator {
 public:
  /*!
  \brief Applies a gradient operator encoded as a dense matrix to a 1D grid.

  \param[in] grad Gradient operator.
  \param[in] grid 1D grid.
  \param[out] out Resulting 1D grid.
  */
  static void ApplyDenseMatrixGradientOn1DGrid(DenseMatrix &grad,
                                               UniStgGrid1D &grid,
                                               UniStgGrid1D &out);

  /*!
  \brief Applies a divergence operator encoded as a dense matrix to a 1D grid.

  \param[in] div Divergence operator.
  \param[in] grid 1D grid.
  \param[out] out Resulting 1D grid.
  */
  static void ApplyDenseMatrixDivergenceOn1DGrid(DenseMatrix &div,
                                                 UniStgGrid1D &grid,
                                                 UniStgGrid1D &out);

  /*!
  \brief Applies a Laplacian operator encoded as a dense matrix to a 1D grid.

  \param[in] lap Laplacian operator.
  \param[in] grid 1D grid.
  \param[out] out Resulting 1D grid.
  */
  static void ApplyDenseMatrixLaplacianOn1DGrid(DenseMatrix &lap,
                                                UniStgGrid1D &grid,
                                                UniStgGrid1D &out);
};
}
#endif  // End of: MTK_INCLUDE_OPERATOR_APPLICATOR_H_
