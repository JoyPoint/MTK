/*!
\file mtk_grad_3d.h

\brief Includes the definition of the class Grad3D.

This class implements a 3D gradient operator, constructed using the
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

#ifndef MTK_INCLUDE_MTK_GRAD_3D_H_
#define MTK_INCLUDE_MTK_GRAD_3D_H_

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_3d.h"

namespace mtk{

/*!
\class Grad3D

\ingroup c07-mim_ops

\brief Implements a 3D mimetic gradient operator.

This class implements a 3D gradient operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).
*/
class Grad3D {
 public:
  /*!
	\brief Default constructor.
	*/
  Grad3D();

  /*!
  \brief Copy constructor.

  \param [in] div Given divergence.
  */
  Grad3D(const Grad3D &grad);

  /*!
  \brief Destructor.
  */
  ~Grad3D();

  /*!
  \brief Factory method implementing the CBS Algorithm to build operator.

  \return Success of the construction.
  */
  bool ConstructGrad3D(const UniStgGrid3D &grid,
                       int order_accuracy = kDefaultOrderAccuracy,
                       Real mimetic_threshold = kDefaultMimeticThreshold);

  /*!
  \brief Return the operator as a dense matrix.

  \return The operator as a dense matrix.
  */
  DenseMatrix ReturnAsDenseMatrix() const;

 private:
  DenseMatrix gradient_;  ///< Actual operator.

  int order_accuracy_;  ///< Order of accuracy.

  Real mimetic_threshold_;  ///< Mimetic Threshold.
};
}
#endif  // End of: MTK_INCLUDE_MTK_GRAD_3D_H_
