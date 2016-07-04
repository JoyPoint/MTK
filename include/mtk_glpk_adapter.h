/*!
\file mtk_glpk_adapter.h

\brief Declaration of an adapter class for the GLPK API.

Declaration of a class that contains a collection of static member functions,
that possess direct access to the underlying structure of the matrices, thus
allowing programmers to exploit some of the numerical methods implemented in
the GLPK.

The **GLPK (GNU Linear Programming Kit)** package is intended for solving
large-scale linear programming (LP), mixed integer programming (MIP), and other
related problems. It is a set of routines written in ANSI C and organized in the
form of a callable library.

\sa http://www.gnu.org/software/glpk/

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

#ifndef MTK_INCLUDE_GLPK_ADAPTER_H_
#define MTK_INCLUDE_GLPK_ADAPTER_H_

#include <iostream>
#include <iomanip>

#include "glpk.h"

#include "mtk_foundations.h"
#include "mtk_dense_matrix.h"

namespace mtk {

/*!
\class GLPKAdapter

\ingroup c05-num_methods

\brief Adapter class for the GLPK API.

This class contains a collection of static member functions, that possess direct
access to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the GLPK.

The **GLPK (GNU Linear Programming Kit)** package is intended for solving
large-scale linear programming (LP), mixed integer programming (MIP), and other
related problems. It is a set of routines written in ANSI C and organized in the
form of a callable library.

\warning We use the GLPK temporarily in order to test the CBSA, but it will be
removed due to potential licensing issues.

\sa http://www.gnu.org/software/glpk/
*/
class GLPKAdapter {
 public:
  /*!
  \brief Solves a CLO problem and compares the solution to a reference solution.

  This routine is the pivot of the CBSA. It solves a Constrained Linear
  Optimization (CLO) problem, and it compares the attained solution to a given
  reference solution. This comparison is done computing the norm-2 relative
  error.

  \param[in] AA Given matrix.
  \param[in] nrows Number of rows of the matrix.
  \param[in] ncols Number of rows of the matrix.
  \param[in] kk Length of the RHS vector of constraints.
  \param[in] hh RHS vector of constraints.
  \param[in,out] qq Output decision vector.
  \param[in] robjective Row of the system to be chosen as objective function.
  \param[in] mimetic_threshold Mimetic threshold.
  \param[in] copy Should we actually copy the results to the output?

  \return Relative error computed between attained solution and provided ref.
  */
  static mtk::Real SolveSimplexAndCompare(mtk::Real *AA,
                                          int nrows,
                                          int ncols,
                                          int kk,
                                          mtk::Real *hh,
                                          mtk::Real *qq,
                                          int robjective,
                                          mtk::Real mimetic_threshold,
                                          int copy) noexcept;
};
}
#endif  // End of: MTK_INCLUDE_MTK_GLPK_ADAPTER_H_
