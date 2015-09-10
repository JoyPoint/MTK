/*!
\file mtk_blas_adapter.h

\brief Adapter class for the BLAS API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the BLAS.

The BLAS (Basic Linear Algebra Subprograms) are routines that provide standard
building blocks for performing basic vector and matrix operations. The Level 1
BLAS perform scalar, vector and vector-vector operations, the Level 2 BLAS
perform matrix-vector operations, and the Level 3 BLAS perform matrix-matrix
operations.

The BLAS can be installed from links given in the See Also section of this page.

\sa http://www.netlib.org/blas/

\sa https://software.intel.com/en-us/non-commercial-software-development

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
*/
/*
Copyright (C) 2015, Computational Science Research Center, San Diego State
University. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Modifications to source code should be reported to: esanchez@mail.sdsu.edu
and a copy of the modified files should be reported once modifications are
completed. Documentation related to said modifications should be included.

2. Redistributions of source code must be done through direct
downloads from the project's GitHub page: http://www.csrc.sdsu.edu/mtk

3. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

4. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

5. Usage of the binary form on proprietary applications shall require explicit
prior written permission from the the copyright holders.

6. Neither the name of the copyright holder nor the names of its contributors
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

#ifndef MTK_INCLUDE_BLAS_ADAPTER_H_
#define	MTK_INCLUDE_BLAS_ADAPTER_H_

#include "mtk_dense_matrix.h"

namespace mtk {

/*!
\class BLASAdapter

\ingroup c05-num_methods

\brief Adapter class for the BLAS API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the BLAS.

The **BLAS (Basic Linear Algebra Subprograms)** are routines that provide
standard building blocks for performing basic vector and matrix operations. The
Level 1 BLAS perform scalar, vector and vector-vector operations, the Level 2
BLAS perform matrix-vector operations, and the Level 3 BLAS perform
matrix-matrix operations.

\sa http://www.netlib.org/blas/
*/
class BLASAdapter {
 public:
  /*!
  \brief Compute the \f$ |\mathbf{x}|_2 \f$ of given array \f$ \mathbf{x} \f$.

  \param[in] in         Input array.
  \param[in] in_length  Length of the array.
  */
  static Real RealNRM2(Real *in, int &in_length);

  /*!
  \brief Real-Arithmetic General (Dense matrices) Matrix-Vector Multiplier.

  Performs

  \f[
    \mathbf{y} := \alpha\mathbf{A}\mathbf{x} + \beta\mathbf{y}
  \f]

  \param[in] alpha  First scalar.
  \param[in] aa     Given matrix.
  \param[in] xx     First vector.
  \param[in] beta   Second scalar.
  \param[in,out] yy Second vector (output).

  \sa http://ejspeiro.github.io/Netlib-and-CPP/
  */
  static void RealDenseMV(Real &alpha,
                          DenseMatrix &aa,
                          Real *xx,
                          Real &beta,
                          Real *yy);

  /*!
  \brief Real-Arithmetic General (Dense matrices) Matrix-Matrix multiplier.

  Performs:

  \f[
  \mathbf{C} := \mathbf{A}\mathbf{B}
  \f]

  \param[in] aa First matrix.
  \param[in] bb Second matrix.

  \sa http://ejspeiro.github.io/Netlib-and-CPP/
  */
  static DenseMatrix RealDenseMM(DenseMatrix &aa, DenseMatrix &bb);
};
}
#endif	// End of: MTK_INCLUDE_BLAS_ADAPTER_H_