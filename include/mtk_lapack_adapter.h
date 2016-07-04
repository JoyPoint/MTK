/*!
\file mtk_lapack_adapter.h

\brief Declaration of an adapter class for the LAPACK API.

Declaration of a class that contains a collection of static member functions,
that possess direct access to the underlying structure of the matrices, thus
allowing programmers to exploit some of the numerical methods implemented in
the LAPACK.

The **LAPACK (Linear Algebra PACKage)** is written in Fortran 90 and provides
routines for solving systems of simultaneous linear equations, least-squares
solutions of linear systems of equations, eigenvalue problems, and singular
value problems.

\sa http://www.netlib.org/lapack/

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

#ifndef MTK_INCLUDE_LAPACK_ADAPTER_H_
#define	MTK_INCLUDE_LAPACK_ADAPTER_H_

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"
#include "mtk_uni_stg_grid_2d.h"

namespace mtk {

/*!
\class LAPACKAdapter

\ingroup c05-num_methods

\brief Adapter class for the LAPACK API

This class contains a collection of static member functions, that possess direct
access to the underlying structure of the matrices, thus allowing programmers to
exploit the numerical methods implemented in the LAPACK.

The **LAPACK (Linear Algebra PACKage)** is written in Fortran 90 and provides
routines for solving systems of simultaneous linear equations, least-squares
solutions of linear systems of equations, eigenvalue problems, and singular
value problems.

\sa http://www.netlib.org/lapack/
*/
class LAPACKAdapter {
 public:
  /*!
  \brief Solves a dense system of linear equations.

  Adapts the MTK to LAPACK's dgesv_ routine.

  \param[in] matrix Input matrix.
  \param[in] rhs    Input right-hand sides vector.

  \exception std::bad_alloc
  */
  static int SolveDenseSystem(mtk::DenseMatrix &mm,
                              mtk::Real *rhs);

  /*!
  \brief Solves a dense system of linear equations.

  Adapts the MTK to LAPACK's dgesv_ routine.

  \param[in] matrix Input matrix.
  \param[in] rr     Input right-hand sides matrix.

  \exception std::bad_alloc
  */
  static int SolveDenseSystem(mtk::DenseMatrix &mm,
                              mtk::DenseMatrix &rr);

  /*!
  \brief Solves a dense system of linear equations.

  Adapts the MTK to LAPACK's dgesv_ routine.

  \param[in] matrix Input matrix.
  \param[in] rhs     Input right-hand side from info on a grid.

  \exception std::bad_alloc
  */
  static int SolveDenseSystem(mtk::DenseMatrix &mm,
                              mtk::UniStgGrid1D &rhs);


  /*!
  \brief Solves a dense system of linear equations.

  Adapts the MTK to LAPACK's dgesv_ routine.

  \param[in] matrix Input matrix.
  \param[in] rhs    Input right-hand side from info on a grid.

  \exception std::bad_alloc
  */
  static int SolveDenseSystem(mtk::DenseMatrix &mm,
                              mtk::UniStgGrid2D &rhs);

  /*!
  \brief Solves overdetermined or underdetermined real linear systems.

  Adapts the MTK to LAPACK's routine.

  \param[in,out] matrix Input matrix.

  \return Success of the solution.

  \exception std::bad_alloc
  */
  static int SolveRectangularDenseSystem(const mtk::DenseMatrix &aa,
                                         mtk::Real *ob_,
                                         int ob_ld_);

  /*!
  \brief Performs a QR factorization on a dense matrix.

  Adapts the MTK to LAPACK's  routine.

  \param[in,out] matrix Input matrix.

  \return Matrix \f$ \mathbf{Q}\f$.

  \exception std::bad_alloc
  */
  static mtk::DenseMatrix QRFactorDenseMatrix(DenseMatrix &matrix);
};
}
#endif  // End of: MTK_INCLUDE_LAPACK_ADAPTER_H_
