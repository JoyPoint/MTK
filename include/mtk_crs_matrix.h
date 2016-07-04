/*!
\file mtk_crs_matrix.h

\brief Declaration of a class for a CRS sparse matrix.

The Compressed Row Storage (CRS) format puts the subsequent nonzeros of the
matrix rows in contiguous memory locations. Assuming we have a nonsymmetric
sparse matrix A, we create 3 vectors: one for floating-point numbers (val), and
the other two for integers (col_ind, row_ptr).

The val vector stores the values of the nonzero elements of the matrix, as they
are traversed in a row-wise fashion. The col_ind vector stores the column
indexes of the elements in the val vector. That is, if \f[$val[k] = a_{ij}$\f]
then \f[$col_ind[k] = j$\f]. The row_ptr vector stores the locations in the val
vector that start a row, that is, if \f[$val[k] = a_{ij}$\f] then \f[$row_ptr[i]
\leq k \leq row_ptr[i + 1]$\f]. By convention, we define \f[$row_ptr[n + 1] =
nnz + 1$\f], where \f[$nnz$\f] is the number of nonzeros in the matrix .

The storage savings for this approach is significant. Instead of storing
\f[$n^2$\f] elements, we need only \f[$2nnz + n + 1$\f] storage locations.

\sa http://netlib.org/linalg/html_templates/node91.html

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\todo Get rid of those pointers to std::vector<> in the member variables.
*/
/*
Copyright (C) 2016, Computational Science Research Center, San Diego State
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

#ifndef MTK_INCLUDE_CRS_MATRIX_
#define MTK_INCLUDE_CRS_MATRIX_

#include <iostream>

#include <vector>

#include "mtk_foundations.h"

namespace mtk {

class CRSMatrix {
public:
  friend std::ostream &operator<<(std::ostream &os, const CRSMatrix &matrix);

  CRSMatrix &operator=(const CRSMatrix &in);
  bool operator==(const CRSMatrix &in);

  CRSMatrix();
  CRSMatrix(const CRSMatrix &mm);
  explicit CRSMatrix(int rank);
  explicit CRSMatrix(int num_rows, int num_cols);
  ~CRSMatrix(void);

  Real GetValue(int row_coord, int col_coord) const;
  CRSMatrix &SetValue(Real val, int row_coord, int col_coord);

  CRSMatrix Add(const CRSMatrix &mm) const;
  std::vector<Real> Multiply(const std::vector<Real> &xx) const;
  CRSMatrix Multiply(const CRSMatrix &mm) const;
  static CRSMatrix Kron(const CRSMatrix &aa, const CRSMatrix &bb);

private:
  void Insert(int index, int row, int col, Real val);
  void Remove(int index, int row);

  std::vector<int> *col_ind_;
  std::vector<int> *row_ptr_;
  std::vector<Real> *val_;

  int num_rows_;
  int num_cols_;
};
}
#endif  // End of: MTK_INCLUDE_CRS_MATRIX_
