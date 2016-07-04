/*!
\file mtk_matrix.cc

\brief Definition of the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.

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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <iostream>

#include "mtk_tools.h"
#include "mtk_matrix.h"

mtk::Matrix::Matrix():
	encoded_operator_(mtk::EncodedOperator::NOOP),
  storage_(mtk::MatrixStorage::DENSE),
  ordering_(mtk::MatrixOrdering::ROW_MAJOR),
  num_rows_(),
  num_cols_(),
  num_values_(),
  leading_dimension_(),
  num_low_diags_(),
  num_upp_diags_(),
  bandwidth_(),
  num_null_(),
  num_non_null_(),
  abs_density_(mtk::kZero),
  abs_sparsity_(mtk::kZero),
  num_zero_(),
  num_non_zero_(),
  rel_density_(mtk::kZero),
  rel_sparsity_(mtk::kZero) {}

mtk::Matrix::Matrix(const Matrix &in):
	encoded_operator_(in.encoded_operator_),
  storage_(in.storage_),
  ordering_(in.ordering_),
  num_rows_(in.num_rows_),
  num_cols_(in.num_cols_),
  num_values_(in.num_values_),
  leading_dimension_(in.leading_dimension_),
  num_low_diags_(in.num_low_diags_),
  num_upp_diags_(in.num_upp_diags_),
  bandwidth_(in.bandwidth_),
  num_null_(in.num_null_),
  num_non_null_(in.num_non_null_),
  abs_density_(in.abs_density_),
  abs_sparsity_(in.abs_sparsity_),
  num_zero_(in.num_zero_),
  num_non_zero_(in.num_non_zero_),
  rel_density_(in.rel_density_),
  rel_sparsity_(in.rel_sparsity_) {}

mtk::Matrix::~Matrix() noexcept {}

mtk::EncodedOperator mtk::Matrix::encoded_operator() const noexcept {

	return encoded_operator_;
}

mtk::MatrixStorage mtk::Matrix::storage() const noexcept {

  return storage_;
}

mtk::MatrixOrdering mtk::Matrix::ordering() const noexcept {

  return ordering_;
}

int mtk::Matrix::num_rows() const noexcept {

  return num_rows_;
}

int mtk::Matrix::num_cols() const noexcept {

  return num_cols_;
}

int mtk::Matrix::num_values() const noexcept {

  return num_values_;
}

int mtk::Matrix::leading_dimension() const noexcept {

  return leading_dimension_;
}

int mtk::Matrix::num_zero() const noexcept {

  return num_zero_;
}

int mtk::Matrix::num_non_zero() const noexcept {

  return num_non_zero_;
}

int mtk::Matrix::num_null() const noexcept {

  return num_null_;
}

int mtk::Matrix::num_non_null() const noexcept {

  return num_non_null_;
}

int mtk::Matrix::num_low_diags() const noexcept {

  return num_low_diags_;
}

int mtk::Matrix::num_upp_diags() const noexcept {

  return num_upp_diags_;
}

int mtk::Matrix::bandwidth() const noexcept {

  return bandwidth_;
}

mtk::Real mtk::Matrix::abs_density() const noexcept {

  return abs_density_;
}

mtk::Real mtk::Matrix::rel_density() const noexcept {

  return rel_density_;
}

mtk::Real mtk::Matrix::abs_sparsity() const noexcept {

  return abs_sparsity_;
}

mtk::Real mtk::Matrix::rel_sparsity() const noexcept {

  return rel_sparsity_;
}

void mtk::Matrix::set_encoded_operator(const mtk::EncodedOperator &in)
	noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	bool aux = (in != mtk::EncodedOperator::NOOP) &&
		(in != mtk::EncodedOperator::GRADIENT) &&
		(in != mtk::EncodedOperator::DIVERGENCE) &&
		(in != mtk::EncodedOperator::INTERPOLATION) &&
		(in != mtk::EncodedOperator::CURL) &&
		(in != mtk::EncodedOperator::LAPLACIAN);

	mtk::Tools::Prevent(aux, __FILE__, __LINE__, __func__);
	#endif

	encoded_operator_ = in;
}

void mtk::Matrix::set_storage(const mtk::MatrixStorage &ss) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(!(ss == mtk::MatrixStorage::DENSE ||
                        ss == mtk::MatrixStorage::BANDED ||
                        ss == mtk::MatrixStorage::CRS),
                      __FILE__, __LINE__, __func__);
  #endif

  storage_ = ss;
}

void mtk::Matrix::set_ordering(const mtk::MatrixOrdering &oo) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  bool aux{oo == mtk::MatrixOrdering::ROW_MAJOR ||
           oo == mtk::MatrixOrdering::COL_MAJOR};
  mtk::Tools::Prevent(!aux, __FILE__, __LINE__, __func__);
  #endif

  ordering_ = oo;

  leading_dimension_ = ComputeLeadingDimension(num_rows_, num_cols_);
}

void mtk::Matrix::set_num_rows(const int &num_rows) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_rows < 1, __FILE__, __LINE__, __func__);
  #endif

  num_rows_ = num_rows;
  num_values_ = ComputeNumValues(num_rows_, num_cols_);
  leading_dimension_ = ComputeLeadingDimension(num_rows_, num_cols_);

  num_null_ = num_values_;
  num_non_null_ = 0;
  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);

  num_zero_ = 0;
  num_non_zero_ = num_values_;
  rel_density_ = ComputeRelDensity(num_non_zero_, num_values_);
  rel_sparsity_ = ComputeRelSparsity(rel_density_);
}

void mtk::Matrix::set_num_cols(const int &num_cols) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_cols < 1, __FILE__, __LINE__, __func__);
  #endif

  num_cols_ = num_cols;
  num_values_ = ComputeNumValues(num_rows_, num_cols_);
  leading_dimension_ = ComputeLeadingDimension(num_rows_, num_cols_);

  num_null_ = num_values_;
  num_non_null_ = 0;
  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);

  num_zero_ = 0;
  num_non_zero_ = num_values_;
  rel_density_ = ComputeRelDensity(num_non_zero_, num_values_);
  rel_sparsity_ = ComputeRelSparsity(rel_density_);
}

void mtk::Matrix::set_num_low_diags(const int &num_low_diags) noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_low_diags < 0, __FILE__, __LINE__, __func__);
	#endif

	num_low_diags_ = num_low_diags;
	bandwidth_ = ComputeBandwidth(num_low_diags_, num_upp_diags_);
}

void mtk::Matrix::set_num_upp_diags(const int &num_upp_diags) noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_upp_diags < 0, __FILE__, __LINE__, __func__);
	#endif

	num_upp_diags_ = num_upp_diags;
	bandwidth_ = ComputeBandwidth(num_low_diags_, num_upp_diags_);
}

void mtk::Matrix::set_num_null(const int &num_null) noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_null < 0, __FILE__, __LINE__, __func__);
	#endif

	num_null_ = num_null;
  num_non_null_ = ComputeNumNonNull(num_values_, num_null_);

  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);;
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);
}

void mtk::Matrix::set_num_zero(const int &num_zero) noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_zero < 0, __FILE__, __LINE__, __func__);
	#endif

  num_zero_ = num_zero;
  num_non_zero_ = ComputeNumNonZero(num_values_, num_zero_);
  rel_density_ = ComputeRelDensity(num_non_zero_, num_values_);
  rel_sparsity_ = ComputeRelSparsity(rel_density_);

  num_null_ = num_null_ - num_zero_;
  num_non_null_ = ComputeNumNonNull(num_values_, num_null_);
  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);
}

void mtk::Matrix::IncreaseNumNull() noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_null_ == num_values_, __FILE__, __LINE__, __func__);
	#endif

  num_null_++;
  num_non_null_ = ComputeNumNonNull(num_values_, num_null_);

  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);;
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);
}

void mtk::Matrix::DecreaseNumNull() noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_null_ == 0, __FILE__, __LINE__, __func__);
  #endif

  num_null_--;
  num_non_null_ = ComputeNumNonNull(num_values_, num_null_);

  abs_density_ = ComputeAbsDensity(num_non_null_, num_values_);
  abs_sparsity_ = ComputeAbsSparsity(abs_density_);
}

void mtk::Matrix::IncreaseNumZero() noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_zero_ == num_values_, __FILE__, __LINE__, __func__);
	#endif

  num_zero_++;
  num_non_zero_ = ComputeNumNonZero(num_values_, num_zero_);

  rel_density_ = ComputeRelDensity(num_non_zero_, num_values_);
  rel_sparsity_ = ComputeRelSparsity(rel_density_);

  DecreaseNumNull();
}

void mtk::Matrix::DecreaseNumZero() noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_zero_ == 0, __FILE__, __LINE__, __func__);
  #endif

  num_zero_--;
  num_non_zero_ = ComputeNumNonZero(num_values_, num_zero_);

  rel_density_ = ComputeRelDensity(num_non_zero_, num_values_);
  rel_sparsity_ = ComputeRelSparsity(rel_density_);

  IncreaseNumNull();
}

int mtk::Matrix::ComputeNumValues(const int &num_rows,
															    const int &num_cols) const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_rows < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cols < 0, __FILE__, __LINE__, __func__);
  #endif

	return num_rows*num_cols;
}

int mtk::Matrix::ComputeLeadingDimension(const int &num_rows,
																				 const int &num_cols) const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(num_rows < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cols < 0, __FILE__, __LINE__, __func__);
  #endif

	return (ordering_ == mtk::MatrixOrdering::ROW_MAJOR)?
	    std::max(1,num_cols): std::max(1,num_rows);
}

int mtk::Matrix::ComputeBandwidth(const int &num_low_diags,
									        	      const int &num_upp_diags) const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_low_diags < 0, __FILE__, __LINE__, __func__);
	mtk::Tools::Prevent(num_upp_diags < 0, __FILE__, __LINE__, __func__);
	#endif

	return num_low_diags + 1 + num_upp_diags;
}

int mtk::Matrix::ComputeNumNonNull(const int &num_values, const int &num_null)
	const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_values < 0, __FILE__, __LINE__, __func__);
	mtk::Tools::Prevent(num_null < 0, __FILE__, __LINE__, __func__);
	#endif

	return num_values - num_null;
}

mtk::Real mtk::Matrix::ComputeAbsDensity(const int &num_non_null,
		                                     const int &num_values) const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_non_null < 0, __FILE__, __LINE__, __func__);
	mtk::Tools::Prevent(num_values < 0, __FILE__, __LINE__, __func__);
	#endif

	return (num_values == mtk::kZero)?
		mtk::kZero: (mtk::Real) num_non_null/num_values;
}

mtk::Real mtk::Matrix::ComputeAbsSparsity(const mtk::Real &absolute_density)
	const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(absolute_density < 0, __FILE__, __LINE__, __func__);
	#endif

	return mtk::kOne - absolute_density;
}

int mtk::Matrix::ComputeNumNonZero(const int &num_values, const int &num_zero)
	const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_values < 0, __FILE__, __LINE__, __func__);
	mtk::Tools::Prevent(num_zero < 0, __FILE__, __LINE__, __func__);
	#endif

	return num_values - num_zero;
}

mtk::Real mtk::Matrix::ComputeRelDensity(const int &num_non_zero,
		                                     const int &num_values) const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(num_non_zero < 0, __FILE__, __LINE__, __func__);
	mtk::Tools::Prevent(num_values < 0, __FILE__, __LINE__, __func__);
	#endif

	return (num_values == mtk::kZero)?
		mtk::kZero: mtk::kOne - (mtk::Real) num_non_zero/num_values;
}

mtk::Real mtk::Matrix::ComputeRelSparsity(const mtk::Real &relative_density)
	const noexcept {

	#ifdef MTK_PERFORM_PREVENTIONS
	mtk::Tools::Prevent(relative_density < 0, __FILE__, __LINE__, __func__);
	#endif

	return mtk::kOne - relative_density;
}
