/*!
\file mtk_bc_desc_2d.h

\brief Enforces boundary conditions in either the operator or the grid.

This class presents an interface for the user to specify boundary conditions
on 2D mimetic operators and the grids they are acting on.

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

#ifndef MTK_INCLUDE_BC_DESC_2D_H_
#define MTK_INCLUDE_BC_DESC_2D_H_

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_2d.h"

namespace mtk{

class BCDesc2D {
 public:
  /*!
  \brief Enforces the condition on the operator represented as matrix.

  \param[in,out] matrix Input operator.
  \param[in] west Pointer to function returning west coefficient at (ii,jj).
  \param[in] east Pointer to function returning east coefficient at (ii,jj).
  \param[in] south Pointer to function returning south coefficient at (ii,jj).
  \param[in] north Pointer to function returning north coefficient at (ii,jj).
  */
  static void ImposeOnOperatorMatrix(const DenseMatrix &matrix,
                                     Real (*west)(int ii, int jj),
                                     Real (*east)(int ii, int jj),
                                     Real (*south)(int ii, int jj),
                                     Real (*north)(int ii, int jj));

/*!
  \brief Enforces the condition on the grid.

  \param[in,out] matrix Input operator.
  \param[in] west Pointer to function returning west coefficient at (xx,yy).
  \param[in] east Pointer to function returning east coefficient at (xx,yy).
  \param[in] south Pointer to function returning south coefficient at (xx,yy).
  \param[in] north Pointer to function returning north coefficient at (xx,yy).
  */
  static void ImposeOnGrid(const UniStgGrid2D &grid,
                           Real (*west)(Real xx, Real yy),
                           Real (*east)(Real xx, Real yy),
                           Real (*south)(Real xx, Real yy),
                           Real (*north)(Real xx, Real yy));
};
}
#endif  // End of: MTK_INCLUDE_BC_DESC_2D_H_
