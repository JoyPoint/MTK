/*!
\file mtk_crs_matrix.cc

\brief Definition of a CRS sparse matrix.

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

#include <iostream>
#include <iomanip>

#include <vector>

#include "mtk_crs_matrix.h"

namespace mtk {

std::ostream &operator<<(std::ostream &stream, const mtk::CRSMatrix &matrix) {

	for (int ii = 1; ii <= matrix.num_rows_; ++ii) {
		for (int jj = 1; jj <= matrix.num_cols_; ++jj) {
			if (jj != 1) {
				stream << std::setw(9);
			}
			stream << matrix.GetValue(ii, jj);
		}
		if (ii < matrix.num_rows_) {
			stream << std::endl;
		}
	}
	return stream;
}
}

bool mtk::CRSMatrix::operator==(const mtk::CRSMatrix &right) {

	bool equality_val = ((val_ == nullptr && right.val_ == nullptr) ||
		(val_ != nullptr && right.val_ != nullptr && *(val_) == *(right.val_)));

	bool equality_col_ind = ((col_ind_ == nullptr && right.col_ind_ == nullptr) ||
		(col_ind_ != nullptr && right.col_ind_ != nullptr &&
			*(col_ind_) == *(right.col_ind_)));

	bool equality_row_ptr = *(row_ptr_) == *(right.row_ptr_);

	return equality_val && equality_col_ind && equality_row_ptr;
}

mtk::CRSMatrix::CRSMatrix():
	col_ind_(),
	row_ptr_(),
	val_(),
	num_rows_(),
	num_cols_() {}

mtk::CRSMatrix::CRSMatrix(const CRSMatrix &in) {}

mtk::CRSMatrix::CRSMatrix(int rank) {

	if (rank < 1 || rank < 1) {
		throw "Matrix dimensions cannot be zero or negative.";
	}

	num_rows_ = rank;
	num_cols_ = rank;

	val_ = nullptr;
	col_ind_ = nullptr;
	row_ptr_ = new std::vector<int>(num_rows_ + 1, 1);
}

mtk::CRSMatrix::CRSMatrix(int num_rows, int num_cols) {

	if (num_rows < 1 || num_cols < 1) {
		throw "Matrix dimensions cannot be zero or negative.";
	}

	num_rows_ = num_rows;
	num_cols_ = num_cols;

	val_ = nullptr;
	col_ind_ = nullptr;
	row_ptr_ = new std::vector<int>(num_rows_ + 1, 1);
}

mtk::CRSMatrix::~CRSMatrix(void) {

	if (val_ != nullptr) {
		delete val_;
	}
	if (col_ind_ != nullptr) {
		delete col_ind_;
	}
	if (row_ptr_ != nullptr) {
		delete row_ptr_;
	}
}

mtk::Real mtk::CRSMatrix::GetValue(int row, int col) const {

	if (row < 1 || col < 1 || row > num_rows_ || col > num_cols_) {
		throw "Coordinates out of range.";
	}

	int actual;

	for (int ii = row_ptr_->at(row - 1) - 1;
		   ii < row_ptr_->at(row) - 1; ++ii) {
		actual = col_ind_->at(ii);
		if (actual == col) {
			return val_->at(ii);
		} else if (actual > col) {
			break;
		}
	}
	return 0.0;
}

mtk::CRSMatrix &mtk::CRSMatrix::SetValue(mtk::Real val, int row, int col) {

	if (row < 1 || col < 1 || row > num_rows_ || col > num_cols_) {
		throw "Coordinations out of range.";
	}

	int pos = row_ptr_->at(row - 1) - 1;
	int actual = -1;

	for (; pos < row_ptr_->at(row) - 1; ++pos) {
		actual = col_ind_->at(pos);

		if (actual == col) {
			break;
		} else if (actual > col) {
			break;
		}
	}

	if (actual != col) {
		if (val != 0) {
			Insert(pos, row, col, val);
		}
	} else if (val == 0) {
		Remove(pos, row);
	} else {
		val_->at(pos) = val;
	}
	return *this;
}

std::vector<mtk::Real> mtk::CRSMatrix::Multiply(const std::vector<mtk::Real> &xx) const {

	if (num_cols_ != static_cast<int>(xx.size())) {
		throw "Cannot multiply: Matrix column count and vector size don't match.";
	}

	std::vector<mtk::Real> result(num_rows_, 0.0);

	for (int ii = 1; ii <= num_rows_; ++ii) {
		for (int jj = 1; jj <= num_cols_; ++jj) {
			result[ii - 1] += GetValue(ii, jj)*xx[jj - 1];
		}
	}
	return result;
}

mtk::CRSMatrix mtk::CRSMatrix::Multiply(const mtk::CRSMatrix &mm) const {

	if (num_cols_ != mm.num_rows_) {
		throw "Left matrix column count and right matrix row count don't match.";
	}

	mtk::CRSMatrix result(num_rows_, mm.num_cols_);

	mtk::Real aa;

	for (int ii = 1; ii <= num_rows_; ++ii) {
		for (int jj = 1; jj <= mm.num_cols_; ++jj) {
			aa = 0.0;
			for (int kk = 1; kk <= num_cols_; ++kk) {
				aa += GetValue(ii, kk)*mm.GetValue(kk, jj);
			}
			result.SetValue(aa, ii, jj);
		}
	}
	return result;
}

mtk::CRSMatrix mtk::CRSMatrix::Kron(const mtk::CRSMatrix &aa,
	                  								const mtk::CRSMatrix &bb) {

	/// \todo Implement Kron using the BLAS.

  // Auxiliary variables:
  auto aux1 = aa.num_rows_*bb.num_rows_;
  auto aux2 = aa.num_cols_*bb.num_cols_;

  mtk::CRSMatrix oo(aux1,aux2); // Output matrix.

  auto mm = aa.num_rows_; // Rows of aa.
  auto nn = aa.num_cols_; // Cols of aa.
  auto pp = bb.num_rows_; // Rows of bb.
  auto qq = bb.num_cols_; // Cols of bb.

  // For each element in the aa matrix...
  for (int ii = 0; ii < mm; ++ii) {
  	for (int jj = 0; jj < nn; ++jj) {
  		mtk::Real aa_factor = aa.GetValue(ii + 1, jj + 1);
  		// ... provided said element is different than 0...
  		if (aa_factor != 0.0) {
  			int oo_row_offset = ii*pp;
  			int oo_col_offset = jj*qq;
  			// ... multiply the entire bb matrix times that element.
  			for (int kk = 0; kk < pp; ++kk) {
					for (int ll = 0; ll < qq; ++ll) {
						mtk::Real bb_factor = bb.GetValue(kk + 1, ll + 1);
						int oo_row = oo_row_offset + kk;
						int oo_col = oo_col_offset + ll;
						oo.SetValue(aa_factor*bb_factor, oo_row + 1, oo_col + 1);
					}
  			}
  		}
  	}
  }

  // oo.matrix_properties_.set_storage(mtk::MatrixStorage::CRS);
  // oo.matrix_properties_.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);

  return oo;
}

mtk::CRSMatrix mtk::CRSMatrix::Add(const mtk::CRSMatrix &mm) const {

	if (num_rows_ != mm.num_rows_ || num_cols_ != mm.num_cols_) {
		throw "Cannot add: matrices dimensions don't match.";
	}

	mtk::CRSMatrix result(num_rows_, num_cols_);

	for (int ii = 1; ii <= num_rows_; ++ii) {
		for (int jj = 1; jj <= num_cols_; ++jj) {
			result.SetValue(GetValue(ii, jj) + mm.GetValue(ii, jj), ii, jj);
		}
	}
	return result;
}

void mtk::CRSMatrix::Insert(int index, int row, int col, mtk::Real val) {

	if (val_ == nullptr) {
		val_ = new std::vector<mtk::Real>(1, val);
		col_ind_ = new std::vector<int>(1, col);
	} else {
		val_->insert(val_->begin() + index, val);
		col_ind_->insert(col_ind_->begin() + index, col);
	}
	for (int ii = row; ii <= num_rows_; ++ii) {
		row_ptr_->at(ii) = row_ptr_->at(ii) + 1;
	}
}

void mtk::CRSMatrix::Remove(int index, int row) {

	val_->erase(val_->begin() + index);
	col_ind_->erase(col_ind_->begin() + index);

	for (int ii = row; ii <= num_rows_; ++ii) {
		row_ptr_->at(ii) = row_ptr_->at(ii) - 1;
	}
}
