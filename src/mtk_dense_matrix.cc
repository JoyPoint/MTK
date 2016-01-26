/*!
\file mtk_dense_matrix.nn

\brief Implements a common dense matrix, using a 1D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.

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

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <typeinfo>

#include <vector>

#include <algorithm>

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"
#include "mtk_tools.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::DenseMatrix &in) {

  int mm{in.matrix_properties_.num_rows()};  // Auxiliary.
  int nn{in.matrix_properties_.num_cols()};  // Auxiliary.
  int output_precision{4};
  int output_width{10};

  if (in.matrix_properties_.ordering() == mtk::MatrixOrdering::COL_MAJOR) {
    std::swap(mm, nn);
  }
  for (int ii = 0; ii < mm; ii++) {
    int offset{ii*nn};
    for (int jj = 0; jj < nn; jj++) {
      mtk::Real value = in.data_[offset + jj];
      stream << std::setprecision(output_precision) <<
        std::setw(output_width) << value;
    }
    stream << std::endl;
  }
  if (in.matrix_properties_.ordering() == mtk::MatrixOrdering::COL_MAJOR) {
    std::swap(mm, nn);
  }
  return stream;
}
}

mtk::DenseMatrix& mtk::DenseMatrix::operator =(const mtk::DenseMatrix &in) {

  if(this == &in) {
    return *this;
  }

  matrix_properties_.set_storage(in.matrix_properties_.storage());

  matrix_properties_.set_ordering(in.matrix_properties_.ordering());

  auto aux = in.matrix_properties_.num_rows();
  matrix_properties_.set_num_rows(aux);

  aux = in.matrix_properties().num_cols();
  matrix_properties_.set_num_cols(aux);

  aux = in.matrix_properties().num_zero();
  matrix_properties_.set_num_zero(aux);

  aux = in.matrix_properties().num_null();
  matrix_properties_.set_num_null(aux);

  auto num_rows = matrix_properties_.num_rows();
  auto num_cols = matrix_properties_.num_cols();

  delete [] data_;

  try {
    data_ = new mtk::Real[num_rows*num_cols];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_, mtk::kZero, sizeof(data_[0])*num_rows*num_cols);

  std::copy(in.data_, in.data_ + num_rows*num_cols, data_);

  return *this;
}

bool mtk::DenseMatrix::operator ==(const DenseMatrix &in) {

  bool ans{true};

  auto mm = in.num_rows();
  auto nn = in.num_cols();

  if (mm != matrix_properties_.num_rows() ||
      nn != matrix_properties_.num_cols()) {
    return false;
  }

  for (int ii = 0; ii < mm && ans; ++ii) {
    for (int jj = 0; jj < nn && ans; ++jj) {
      ans = ans &&
        abs(data_[ii*nn + jj] - in.data()[ii*nn + jj]) < mtk::kDefaultTolerance;
    }
  }
  return ans;
}

mtk::DenseMatrix::DenseMatrix(): data_(nullptr) {

  matrix_properties_.set_storage(mtk::MatrixStorage::DENSE);
  matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
}

mtk::DenseMatrix::DenseMatrix(const mtk::DenseMatrix &in) {

  matrix_properties_.set_storage(in.matrix_properties_.storage());

  matrix_properties_.set_ordering(in.matrix_properties_.ordering());

  auto aux = in.matrix_properties_.num_rows();
  matrix_properties_.set_num_rows(aux);

  aux = in.matrix_properties().num_cols();
  matrix_properties_.set_num_cols(aux);

  aux = in.matrix_properties().num_zero();
  matrix_properties_.set_num_zero(aux);

  aux = in.matrix_properties().num_null();
  matrix_properties_.set_num_null(aux);

  auto num_rows = in.matrix_properties_.num_rows();
  auto num_cols = in.matrix_properties_.num_cols();

  try {
    data_ = new mtk::Real[num_rows*num_cols];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_, mtk::kZero, sizeof(data_[0])*num_rows*num_cols);

  std::copy(in.data_,in.data_ + num_rows*num_cols,data_);
}

mtk::DenseMatrix::DenseMatrix(const int &num_rows, const int &num_cols) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_rows < 1, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cols < 1, __FILE__, __LINE__, __func__);
  #endif

  matrix_properties_.set_storage(mtk::MatrixStorage::DENSE);
  matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
  matrix_properties_.set_num_rows(num_rows);
  matrix_properties_.set_num_cols(num_cols);

  try {
    data_ = new mtk::Real[num_rows*num_cols];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_, mtk::kZero, sizeof(data_[0])*num_rows*num_cols);
}

mtk::DenseMatrix::DenseMatrix(const int &rank,
                              const bool &padded,
                              const bool &transpose) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(rank < 1, __FILE__, __LINE__, __func__);
  #endif

  int aux{};  // Used to control the padding.

  if (padded) {
    aux = 1;
  }

  matrix_properties_.set_storage(mtk::MatrixStorage::DENSE);
  matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
  matrix_properties_.set_num_rows(aux + rank + aux);
  matrix_properties_.set_num_cols(rank);

  try {
    data_ = new mtk::Real[matrix_properties_.num_values()];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_,
         mtk::kZero,
         sizeof(data_[0])*(matrix_properties_.num_values()));

  for (auto ii =0 ; ii < matrix_properties_.num_rows(); ++ii) {
    for (auto jj = 0; jj < matrix_properties_.num_cols(); ++jj) {
      data_[ii*matrix_properties_.num_cols() + jj] =
        (ii == jj + aux)? mtk::kOne: mtk::kZero;
    }
  }
  if (transpose) {
    Transpose();
  }
}

mtk::DenseMatrix::DenseMatrix(const mtk::Real *const gen,
                              const int &gen_length,
                              const int &pro_length,
                              const bool &transpose) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(gen == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(gen_length < 1, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(pro_length < 1, __FILE__, __LINE__, __func__);
  #endif

  matrix_properties_.set_storage(mtk::MatrixStorage::DENSE);
  matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
  if (!transpose) {
    matrix_properties_.set_num_rows(gen_length);
    matrix_properties_.set_num_cols(pro_length);
  } else {
    matrix_properties_.set_num_rows(pro_length);
    matrix_properties_.set_num_cols(gen_length);
  }

  int mm = matrix_properties_.num_rows(); // Used to construct this matrix.
  int nn = matrix_properties_.num_cols(); // Used to construct this matrix.

  try {
    data_ = new mtk::Real[mm*nn];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_, mtk::kZero, sizeof(data_[0])*mm*nn);

  if (!transpose) {
    for (auto ii = 0; ii < mm; ii++) {
      for (auto jj = 0; jj < nn; jj++) {
        data_[ii*nn + jj] = pow(gen[ii], (double) jj);
      }
    }
  } else {
    for (auto ii = 0; ii < mm; ii++) {
      for (auto jj = 0; jj < nn; jj++) {
        data_[ii*nn + jj] = pow(gen[jj], (double) ii);
      }
    }
  }
}

mtk::DenseMatrix::~DenseMatrix() {

  delete [] data_;
  data_ = nullptr;
}

mtk::Matrix mtk::DenseMatrix::matrix_properties() const noexcept {

  return matrix_properties_;
}

void mtk::DenseMatrix::SetOrdering(mtk::MatrixOrdering oo) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(!(oo == mtk::MatrixOrdering::ROW_MAJOR || oo ==
mtk::MatrixOrdering::COL_MAJOR),
                      __FILE__, __LINE__, __func__);
  #endif

  matrix_properties_.set_ordering(oo);
}

int mtk::DenseMatrix::num_rows() const noexcept {

  return matrix_properties_.num_rows();
}

int mtk::DenseMatrix::num_cols() const noexcept {

  return matrix_properties_.num_cols();
}

mtk::Real* mtk::DenseMatrix::data() const noexcept {

  return data_;
}

mtk::Real mtk::DenseMatrix::GetValue(
    const int &mm,
    const int &nn) const noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(mm < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 0, __FILE__, __LINE__, __func__);
  #endif

  return data_[mm*matrix_properties_.num_cols() + nn];
}

void  mtk::DenseMatrix::SetValue(
    const int &mm,
    const int &nn,
    const mtk::Real &val) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(mm < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 0, __FILE__, __LINE__, __func__);
  #endif

  data_[mm*matrix_properties_.num_cols() + nn] = val;
}

void mtk::DenseMatrix::Transpose() {

  /// \todo Improve this so that no extra arrays have to be created.

  mtk::Real *data_transposed{}; // Buffer.

  int mm = matrix_properties_.num_rows(); // Used to construct this matrix.
  int nn = matrix_properties_.num_cols(); // Used to construct this matrix.

  try {
    data_transposed = new mtk::Real[mm*nn];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(data_transposed,
         mtk::kZero,
         sizeof(data_transposed[0])*mm*nn);

  // Assign the values to their transposed position.
  for (auto ii = 0; ii < mm; ++ii) {
    for (auto jj = 0; jj < nn; ++jj) {
      data_transposed[jj*mm + ii] = data_[ii*nn + jj];
    }
  }

  // Swap pointers.
  auto tmp = data_; // Temporal holder.
  data_ = data_transposed;
  delete [] tmp;
  tmp = nullptr;

  matrix_properties_.set_num_rows(nn);
  matrix_properties_.set_num_cols(mm);
}

void mtk::DenseMatrix::OrderRowMajor() {

  if (matrix_properties_.ordering() == mtk::MatrixOrdering::COL_MAJOR) {

    /// \todo Improve this so that no new arrays have to be created.

    mtk::Real *data_transposed{}; // Buffer.

    int mm = matrix_properties_.num_rows(); // Used to construct this matrix.
    int nn = matrix_properties_.num_cols(); // Used to construct this matrix.

    try {
      data_transposed = new mtk::Real[mm*nn];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(data_transposed,
          mtk::kZero,
          sizeof(data_transposed[0])*mm*nn);

    // Assign the values to their transposed position.
    std::swap(mm, nn);
    for (auto ii = 0; ii < mm; ++ii) {
      for (auto jj = 0; jj < nn; ++jj) {
        data_transposed[jj*mm + ii] = data_[ii*nn + jj];
      }
    }
    std::swap(mm, nn);

    // Swap pointers.
    auto tmp = data_; // Temporal holder.
    data_ = data_transposed;
    delete [] tmp;
    tmp = nullptr;

    matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
  }
}

void mtk::DenseMatrix::OrderColMajor() {

  if (matrix_properties_.ordering() == mtk::MatrixOrdering::ROW_MAJOR) {

    /// \todo Improve this so that no new arrays have to be created.

    mtk::Real *data_transposed{}; // Buffer.

    int mm = matrix_properties_.num_rows(); // Used to construct this matrix.
    int nn = matrix_properties_.num_cols(); // Used to construct this matrix.

    try {
      data_transposed = new mtk::Real[mm*nn];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(data_transposed,
          mtk::kZero,
          sizeof(data_transposed[0])*mm*nn);

    // Assign the values to their transposed position.
    for (auto ii = 0; ii < mm; ++ii) {
      for (auto jj = 0; jj < nn; ++jj) {
        data_transposed[jj*mm + ii] = data_[ii*nn + jj];
      }
    }

    // Swap pointers.
    auto tmp = data_; // Temporal holder.
    data_ = data_transposed;
    delete [] tmp;
    tmp = nullptr;

    matrix_properties_.set_ordering(mtk::MatrixOrdering::COL_MAJOR);
  }
}

mtk::DenseMatrix mtk::DenseMatrix::Kron(const mtk::DenseMatrix &aa,
                                        const mtk::DenseMatrix &bb) {

  /// \todo Implement Kron using the BLAS.

  int row_offset{}; // Offset for rows.
  int col_offset{}; // Offset for rows.

  mtk::Real aa_factor{};    // Used in computation.

  // Auxiliary variables:
  auto aux1 = aa.matrix_properties_.num_rows()*bb.matrix_properties_.num_rows();
  auto aux2 = aa.matrix_properties_.num_cols()*bb.matrix_properties_.num_cols();

  mtk::DenseMatrix output(aux1,aux2); // Output matrix.

  int kk_num_cols{output.matrix_properties_.num_cols()}; // Aux.

  auto mm = aa.matrix_properties_.num_rows(); // Rows of aa.
  auto nn = aa.matrix_properties_.num_cols(); // Cols of aa.
  auto pp = bb.matrix_properties_.num_rows(); // Rows of bb.
  auto qq = bb.matrix_properties_.num_cols(); // Cols of bb.

  for (auto ii = 0; ii < mm; ++ii) {
    row_offset = ii*pp;
    for (auto jj = 0; jj < nn; ++jj) {
      col_offset = jj*qq;
      aa_factor = aa.data_[ii*nn + jj];
      for (auto ll = 0; ll < pp; ++ll) {
        for (auto oo = 0; oo < qq; ++oo) {
          auto index = (ll + row_offset)*kk_num_cols + (oo + col_offset);
          output.data_[index] = aa_factor*bb.data_[ll*qq + oo];
        }
      }
    }
  }

  output.matrix_properties_.set_storage(mtk::MatrixStorage::DENSE);
  output.matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);

  return output;
}

bool mtk::DenseMatrix::WriteToFile(const std::string &filename) const {

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open(filename);

  if (!output_dat_file.is_open()) {
    return false;
  }

  int mm{matrix_properties_.num_rows()};
  int nn{matrix_properties_.num_cols()};

  for (int ii = 0; ii < mm; ++ii) {
    int offset{ii*nn};
    for (int jj = 0; jj < nn; ++jj) {
      output_dat_file << ii << ' ' << jj << ' ' << data_[offset + jj] <<
        std::endl;
    }
  }

  output_dat_file.close();

  return true;
}
