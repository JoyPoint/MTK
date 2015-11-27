/*!
\file mtk_matrix.cc

\brief Implementing the representation of a matrix in the MTK.

Implementation of the representation for the matrices implemented in the MTK.

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
#include <cstdio>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <iostream>

#include "mtk_tools.h"
#include "mtk_matrix.h"

mtk::Matrix::Matrix():
  storage_(mtk::DENSE),
  ordering_(mtk::ROW_MAJOR),
  num_rows_(),
  num_cols_(),
  num_values_(),
  ld_(),
  num_zero_(),
  num_non_zero_(),
  num_null_(),
  num_non_null_(),
  kl_(),
  ku_(),
  bandwidth_(),
  abs_density_(),
  rel_density_(),
  abs_sparsity_(),
  rel_sparsity_() {}

mtk::Matrix::Matrix(const Matrix &in):
  storage_(in.storage_),
  ordering_(in.ordering_),
  num_rows_(in.num_rows_),
  num_cols_(in.num_cols_),
  num_values_(in.num_values_),
  ld_(in.ld_),
  num_zero_(in.num_zero_),
  num_non_zero_(in.num_non_zero_),
  num_null_(in.num_null_),
  num_non_null_(in.num_non_null_),
  kl_(in.kl_),
  ku_(in.ku_),
  bandwidth_(in.bandwidth_),
  abs_density_(in.abs_density_),
  rel_density_(in.rel_density_),
  abs_sparsity_(in.abs_sparsity_),
  rel_sparsity_(in.rel_sparsity_) {}

mtk::Matrix::~Matrix() noexcept {}

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

int mtk::Matrix::ld() const noexcept {

  return ld_;
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

int mtk::Matrix::kl() const noexcept {

  return kl_;
}

int mtk::Matrix::ku() const noexcept {

  return ku_;
}

int mtk::Matrix::bandwidth() const noexcept {

  return bandwidth_;
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

void mtk::Matrix::set_storage(const mtk::MatrixStorage &ss) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(!(ss == mtk::DENSE ||
                        ss == mtk::BANDED ||
                        ss == mtk::CRS),
                      __FILE__, __LINE__, __func__);
  #endif

  storage_ = ss;
}

void mtk::Matrix::set_ordering(const mtk::MatrixOrdering &oo) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(!(oo == mtk::ROW_MAJOR || oo == mtk::COL_MAJOR),
                      __FILE__, __LINE__, __func__);
  #endif

  ordering_ = oo;

  ld_ = (ordering_ == mtk::ROW_MAJOR)?
    std::max(1,num_cols_): std::max(1,num_rows_);
}

void mtk::Matrix::set_num_rows(const int &in) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(in < 1, __FILE__, __LINE__, __func__);
  #endif

  num_rows_ = in;
  num_values_ = num_rows_*num_cols_;
  ld_ = (ordering_ == mtk::ROW_MAJOR)?
    std::max(1,num_cols_): std::max(1,num_rows_);
}

void mtk::Matrix::set_num_cols(const int &in) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(in < 1, __FILE__, __LINE__, __func__);
  #endif

  num_cols_ = in;
  num_values_ = num_rows_*num_cols_;
  ld_ = (ordering_ == mtk::ROW_MAJOR)?
    std::max(1,num_cols_): std::max(1,num_rows_);
}

void mtk::Matrix::set_num_zero(const int &in) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(in < 0, __FILE__, __LINE__, __func__);
  #endif

  num_zero_ = in;
  num_non_zero_ = num_values_ - num_zero_;

  /// \bug -nan assigned on construction time due to num_values_ being 0.
  rel_density_ = (mtk::Real) num_non_zero_/num_values_;
  rel_sparsity_ = 1.0 - rel_density_;
}

void mtk::Matrix::set_num_null(const int &in) noexcept {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(in < 0, __FILE__, __LINE__, __func__);
  #endif

  num_null_ = in;
  num_non_null_ = num_values_ - num_null_;

  /// \bug -nan assigned on construction time due to num_values_ being 0.
  abs_density_ = (mtk::Real) num_non_null_/num_values_;
  abs_sparsity_ = 1.0 - abs_density_;
}

void mtk::Matrix::IncreaseNumZero() noexcept {

  /// \todo Review the definition of sparse matrices properties.

  num_zero_++;
  num_non_zero_ = num_values_ - num_zero_;
  rel_density_ = (mtk::Real) num_non_zero_/num_values_;
  rel_sparsity_ = 1.0 - rel_density_;
}

void mtk::Matrix::IncreaseNumNull() noexcept {

  /// \todo Review the definition of sparse matrices properties.

  num_null_++;
  num_non_null_ = num_values_ - num_null_;
  abs_density_ = (mtk::Real) num_non_null_/num_values_;
  abs_sparsity_ = 1.0 - abs_density_;
}
