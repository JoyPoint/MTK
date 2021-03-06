/*!
\file mtk_grad_1d.cc

\brief Definition of the class Grad1D.

This class implements a 1D gradient matrix operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\todo Overload ostream operator as in mtk::Lap1D.

\todo Implement creation of \f$ \mathbf{\Lambda}\f$ w. mtk::BLASAdapter.
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

#include <cmath>
#include <cstring>

#include <iostream>
#include <iomanip>

#ifdef MTK_VERBOSE_WEIGHTS
#include <fstream>
#endif

#include <limits>
#include <algorithm>

#include "mtk_tools.h"
#include "mtk_blas_adapter.h"
#include "mtk_lapack_adapter.h"
#include "mtk_glpk_adapter.h"
#include "mtk_grad_1d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::Grad1D &in) {

  int output_precision{4};
  int output_width{8};

  /// 1. Print order of accuracy.

  stream << "Order of accuracy: " << in.gradient_[0] << std::endl;

  /// 2. Print approximating coefficients for the interior.

  stream << "Interior stencil: " << std::endl;
  for (auto ii = 1; ii <= in.order_accuracy_; ++ii) {
    stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.gradient_[ii] << ' ';
  }
  stream << std::endl;

  /// 3. Print mimetic weights.

  stream << "Weights:" << std::endl;
  for (auto ii = in.order_accuracy_ + 1; ii <= 2*in.order_accuracy_; ++ii) {
    stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.gradient_[ii] << ' ';
  }
  stream << std::endl;

  /// 4. Print mimetic approximations at the boundary.

  int offset{2*in.order_accuracy_ + 1};
  int mm {};
  if (in.order_accuracy_ > mtk::kDefaultOrderAccuracy) {
    for (auto ii = 0; ii < in.num_bndy_approxs_ ; ++ii) {
      stream << "Boundary row " << ii + 1 << ":" << std::endl;
      for (auto jj = 0; jj < in.num_bndy_coeffs_; jj++) {
        auto value = in.gradient_[offset + (mm)];
        stream << std::setprecision(output_precision) <<
        std::setw(output_width) << value << ' ';
        mm++;
      }
      stream << std::endl;
      stream << "Sum of elements in boundary row " << ii + 1 << ": " <<
        in.sums_rows_mim_bndy_[ii];
      stream << std::endl;
    }
  } else {
    stream << "Boundary row 1:" << std::endl;
    stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.gradient_[offset + 0] << ' ';
    stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.gradient_[offset + 1] << ' ';
    stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.gradient_[offset + 2] << ' ';
    stream << std::endl;
    stream << "Sum of elements in boundary row 1: " <<
      in.gradient_[offset + 0] + in.gradient_[offset + 1] +
        in.gradient_[offset + 2];
    stream << std::endl;
  }

  return stream;
}
}

mtk::Grad1D::Grad1D():
  order_accuracy_(mtk::kDefaultOrderAccuracy),
  dim_null_(),
  num_bndy_approxs_(),
  num_bndy_coeffs_(),
  gradient_length_(),
  minrow_(),
  row_(),
  num_feasible_sols_(),
  coeffs_interior_(),
  prem_apps_(),
  weights_crs_(),
  weights_cbs_(),
  mim_bndy_(),
  gradient_(),
  mimetic_threshold_(mtk::kDefaultMimeticThreshold),
  mimetic_measure_(mtk::kZero),
  sums_rows_mim_bndy_() {}

mtk::Grad1D::Grad1D(const Grad1D &grad):
  order_accuracy_(grad.order_accuracy_),
  dim_null_(grad.dim_null_),
  num_bndy_approxs_(grad.num_bndy_approxs_),
  num_bndy_coeffs_(grad.num_bndy_coeffs_),
  gradient_length_(grad.gradient_length_),
  minrow_(grad.minrow_),
  row_(grad.row_),
  num_feasible_sols_(grad.num_feasible_sols_),
  coeffs_interior_(grad.coeffs_interior_),
  prem_apps_(grad.prem_apps_),
  weights_crs_(grad.weights_crs_),
  weights_cbs_(grad.weights_cbs_),
  mim_bndy_(grad.mim_bndy_),
  gradient_(grad.gradient_),
  mimetic_threshold_(grad.mimetic_threshold_),
  mimetic_measure_(grad.mimetic_measure_),
  sums_rows_mim_bndy_(grad.sums_rows_mim_bndy_) {}

mtk::Grad1D::~Grad1D() {

  delete[] coeffs_interior_;
  coeffs_interior_ = nullptr;

  delete[] prem_apps_;
  prem_apps_ = nullptr;

  delete[] weights_crs_;
  weights_crs_ = nullptr;

  delete[] weights_cbs_;
  weights_cbs_ = nullptr;

  delete[] mim_bndy_;
  mim_bndy_ = nullptr;

  delete[] gradient_;
  gradient_ = nullptr;
}

bool mtk::Grad1D::ConstructGrad1D(int order_accuracy, Real mimetic_threshold) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(order_accuracy < 2, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent((order_accuracy%2) != 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(mimetic_threshold <= mtk::kZero,
                      __FILE__, __LINE__, __func__);

  if (order_accuracy >= mtk::kCriticalOrderAccuracyGrad) {
    std::cout << "WARNING: Numerical accuracy is high." << std::endl;
  }

  std::cout << "order_accuracy_ = " << order_accuracy << std::endl;
  std::cout << "mimetic_threshold_ = " << mimetic_threshold << std::endl;
  #endif

  order_accuracy_ = order_accuracy;
  mimetic_threshold_ = mimetic_threshold;

  /// 1. Compute stencil for the interior cells.
  bool abort_construction = ComputeStencilInteriorGrid();

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!abort_construction) {
    std::cerr << "Could NOT complete stage 1." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  #endif

  // At this point, we already have the values for the interior stencil stored
  // in the coeffs_interior_ array.

  dim_null_ = order_accuracy_/2 - 1;

  num_bndy_approxs_ = dim_null_ + 1;

  #ifdef MTK_PRECISION_DOUBLE
  num_bndy_coeffs_ = (int) (3.0*((mtk::Real) order_accuracy_)/2.0);
  #else
  num_bndy_coeffs_ = (int) (3.0f*((mtk::Real) order_accuracy_)/2.0f);
  #endif

  /// 2. Compute a rational null-space from the first matrix transposed.

  // For this we will follow recommendations given in:
  //
  // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=4506
  //
  // We will compute the QR Factorization of the transpose, as in the
  // following (MATLAB) pseudo-code:
  //
  // [Q,R] = qr(V'); % Full QR as defined in
  // % http://www.stanford.edu/class/ee263/notes/qr_matlab.pdf
  //
  // null-space = Q(:, last (order_accuracy_/2 - 1) columns of Q );
  //
  // However, given the nature of the Vandermonde matrices we've just
  // computed, they all posses the same null-space. Therefore, we impose the
  // convention of computing the null-space of the first Vandermonde matrix
  // (west boundary).

  // In the case of the gradient, the first Vandermonde system has a unique
  // solution for the case of second-order-accuracy. Ergo, the Vandermonde
  // matrix used to assemble said system, will have an empty null-space.

  // Therefore, we only compute a rational basis for the case of order higher
  // than second.

  if (dim_null_ > 0) {

    abort_construction = ComputeRationalBasisNullSpace();

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!abort_construction) {
      std::cerr << "Could NOT complete stage 2.1." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      return false;
    }
    #endif
  }

  /// 3. Compute preliminary approximation (non-mimetic) on the boundaries.
  abort_construction = ComputePreliminaryApproximations();

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!abort_construction) {
    std::cerr << "Could NOT complete stage 2.2." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  #endif

  /// 4. Compute quadrature weights to impose the mimetic conditions.
  abort_construction = ComputeWeights();

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!abort_construction) {
    std::cerr << "Could NOT complete stage 2.3." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  #endif

  /// 5. Compute real approximation (mimetic) on the boundaries.
  if (dim_null_ > 0) {

    abort_construction = ComputeStencilBoundaryGrid();

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!abort_construction) {
      std::cerr << "Could NOT complete stage 2.4." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      return false;
    }
    #endif
  }

  /// 6. Assemble operator.

  // Once we have the following three collections of data:
  //   (a) the coefficients for the interior,
  //   (b) the coefficients for the boundary (if it applies),
  //   (c) and the weights (if it applies),
  // we will store everything in the output array:

  abort_construction = AssembleOperator();

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!abort_construction) {
    std::cerr << "Could NOT complete stage 3." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  #endif

  return true;
}

int mtk::Grad1D::num_bndy_coeffs() const {

  return num_bndy_coeffs_;
}

mtk::Real *mtk::Grad1D::coeffs_interior() const {

  return coeffs_interior_;
}

mtk::Real *mtk::Grad1D::weights_crs() const {

  return weights_crs_;
}

mtk::Real *mtk::Grad1D::weights_cbs() const {

  return weights_cbs_;
}

int mtk::Grad1D::num_feasible_sols() const {

  return num_feasible_sols_;
}

mtk::DenseMatrix mtk::Grad1D::mim_bndy() const {

  mtk::DenseMatrix xx(dim_null_ + 1, 3*order_accuracy_/2);

  auto counter = 0;
  for (auto ii = 0; ii < dim_null_ + 1; ++ii) {
    for(auto jj = 0; jj < 3*order_accuracy_/2; ++jj) {
      xx.SetValue(ii,jj, gradient_[2*order_accuracy_ + 1 + counter]);
      counter++;
    }
  }

  return xx;
}

std::vector<mtk::Real> mtk::Grad1D::sums_rows_mim_bndy() const {

  return sums_rows_mim_bndy_;
}

mtk::Real mtk::Grad1D::mimetic_measure() const {

  return mimetic_measure_;
}

mtk::DenseMatrix mtk::Grad1D::ReturnAsDenseMatrix(mtk::Real west,
                                                  mtk::Real east,
                                                  int num_cells_x) const {

  int nn{num_cells_x}; // Number of cells on the grid.

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(east < west, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 3*order_accuracy_ - 2, __FILE__, __LINE__, __func__);
  #endif

  mtk::Real delta_x = (east - west)/((mtk::Real) num_cells_x);

  mtk::Real inv_delta_x{mtk::kOne/delta_x};

  int gg_num_rows = nn + 1;
  int gg_num_cols = nn + 2;
  int num_extra_rows = order_accuracy_/2;
  int elements_per_extra_row = num_bndy_coeffs_;

  // Output matrix featuring sizes for gradient operators.
  mtk::DenseMatrix out(gg_num_rows, gg_num_cols);

  out.set_encoded_operator(mtk::EncodedOperator::GRADIENT);

  /// 1. Insert mimetic boundary at the west.

  auto ee_index = 0;
  for (auto ii = 0; ii < num_extra_rows; ii++) {
    auto cc = 0;
    for(auto jj = 0 ; jj < gg_num_cols; jj++) {
      if(cc >= elements_per_extra_row) {
        out.SetValue(ii, jj, mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                     gradient_[2*order_accuracy_ + 1 + ee_index++]*inv_delta_x);
        cc++;
      }
    }
  }

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = num_extra_rows; ii < gg_num_rows - num_extra_rows; ii++) {
    auto jj = ii - num_extra_rows + 1;
    for (auto cc = 0; cc < order_accuracy_; cc++, jj++) {
      out.SetValue(ii, jj, coeffs_interior_[cc]*inv_delta_x);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the mimetic boundaries.

  ee_index = 0;
  for (auto ii = gg_num_rows - 1; ii >= gg_num_rows - num_extra_rows; ii--) {
    auto cc = 0;
    for (auto jj = gg_num_cols - 1; jj >= 0; jj--) {
      if(cc >= elements_per_extra_row) {
        out.SetValue(ii,jj,mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                     -gradient_[2*order_accuracy_ + 1 +
ee_index++]*inv_delta_x);
        cc++;
      }
     }
  }

  return out;
}

mtk::DenseMatrix mtk::Grad1D::ReturnAsDenseMatrix(
  const UniStgGrid1D &grid) const {

  int nn{grid.num_cells_x()}; // Number of cells on the grid.

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 3*order_accuracy_ - 2, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(grid.field_nature() != mtk::FieldNature::SCALAR,
                      __FILE__, __LINE__, __func__);
  #endif

  mtk::Real inv_delta_x{mtk::kOne/grid.delta_x()};

  int gg_num_rows = nn + 1;
  int gg_num_cols = nn + 2;
  int num_extra_rows = order_accuracy_/2;
  int elements_per_row = num_bndy_coeffs_;

  // Output matrix featuring sizes for gradient operators.
  mtk::DenseMatrix out(gg_num_rows, gg_num_cols);

  out.set_encoded_operator(mtk::EncodedOperator::GRADIENT);

  /// 1. Insert mimetic boundary at the west.

  auto ee_index = 0;
  for (auto ii = 0; ii < num_extra_rows; ii++) {
    auto cc = 0;
    for(auto jj = 0 ; jj < gg_num_cols; jj++) {
      if(cc >= elements_per_row) {
        out.SetValue(ii, jj, mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                     gradient_[2*order_accuracy_ + 1 + ee_index++]*inv_delta_x);
        cc++;
      }
    }
  }

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = num_extra_rows; ii < gg_num_rows - num_extra_rows; ii++) {
    auto jj = ii - num_extra_rows + 1;
    for (auto cc = 0; cc < order_accuracy_; cc++, jj++) {
      out.SetValue(ii, jj, coeffs_interior_[cc]*inv_delta_x);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the mimetic boundaries.

  ee_index = 0;
  for (auto ii = gg_num_rows - 1; ii >= gg_num_rows - num_extra_rows; ii--) {
    auto cc = 0;
    for (auto jj = gg_num_cols - 1; jj >= 0; jj--) {
      if(cc >= elements_per_row) {
        out.SetValue(ii,jj,mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                    -gradient_[2*order_accuracy_ + 1 + ee_index++]*inv_delta_x);
        cc++;
      }
     }
  }

  return out;
}

mtk::DenseMatrix mtk::Grad1D::ReturnAsDimensionlessDenseMatrix(
  int num_cells_x) const {

  int nn{num_cells_x}; // Number of cells on the grid.

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 3*order_accuracy_ - 2, __FILE__, __LINE__, __func__);
  #endif

  int gg_num_rows = nn + 1;
  int gg_num_cols = nn + 2;
  int elements_per_extra_row = num_bndy_coeffs_;
  int num_extra_rows = order_accuracy_/2;

  // Output matrix featuring sizes for gradient operators.
  mtk::DenseMatrix out(gg_num_rows, gg_num_cols);

  out.set_encoded_operator(mtk::EncodedOperator::GRADIENT);

  /// 1. Insert mimetic boundary at the west.

  auto ee_index = 0;
  for (auto ii = 0; ii < num_extra_rows; ii++) {
    auto cc = 0;
    for(auto jj = 0 ; jj < gg_num_cols; jj++) {
      if(cc >= elements_per_extra_row) {
        out.SetValue(ii, jj, mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                     gradient_[2*order_accuracy_ + 1 + ee_index++]);
        cc++;
      }
    }
  }

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = num_extra_rows; ii < gg_num_rows - num_extra_rows; ii++) {
    auto jj = ii - num_extra_rows + 1;
    for (auto cc = 0; cc < order_accuracy_; cc++, jj++) {
      out.SetValue(ii, jj, coeffs_interior_[cc]);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the mimetic boundaries.

  ee_index = 0;
  for (auto ii = gg_num_rows - 1; ii >= gg_num_rows - num_extra_rows; ii--) {
    auto cc = 0;
    for (auto jj = gg_num_cols - 1; jj >= 0; jj--) {
      if(cc >= elements_per_extra_row) {
        out.SetValue(ii,jj,mtk::kZero);
      } else {
        out.SetValue(ii,jj,
                     -gradient_[2*order_accuracy_ + 1 + ee_index++]);
        cc++;
      }
     }
  }

  return out;
}

bool mtk::Grad1D::ComputeStencilInteriorGrid() {

  /// 1. Create vector for interior spatial coordinates.

  mtk::Real* pp{}; // Spatial coordinates to create interior stencil.

  try {
    pp = new mtk::Real[order_accuracy_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(pp, mtk::kZero, sizeof(pp[0])*order_accuracy_);

  #ifdef MTK_PRECISION_DOUBLE
  pp[0] = 1.0/2.0 - ((mtk::Real) order_accuracy_)/2.0;
  #else
  pp[0] = 1.0f/2.0f - ((mtk::Real) order_accuracy_)/2.0f;
  #endif

  for (auto ii = 1; ii < order_accuracy_; ++ii) {
    pp[ii] = pp[ii - 1] + mtk::kOne;
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "pp =" << std::endl;
  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    std::cout << std::setw(12) << pp[ii];
  }
  std::cout << std::endl << std::endl;
  #endif

  /// 2. Create Vandermonde matrix (using interior coordinates as generator).

  bool transpose{false};

  mtk::DenseMatrix vander_matrix(pp,order_accuracy_,order_accuracy_,transpose);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "vander_matrix = " << std::endl;
  std::cout << vander_matrix << std::endl << std::endl;
  #endif

  /// 3. Create order-selector vector.

  try {
    coeffs_interior_ = new mtk::Real[order_accuracy_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(coeffs_interior_, mtk::kZero,
sizeof(coeffs_interior_[0])*order_accuracy_);

  coeffs_interior_[1] = mtk::kOne;

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "oo =" << std::endl;
  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    std::cout << std::setw(12) << coeffs_interior_[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 4. Solve dense Vandermonde system to attain the interior coefficients.

  int info{mtk::LAPACKAdapter::SolveDenseSystem(vander_matrix,
                                                coeffs_interior_)};

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cout << "System solved! Interior stencil attained!" << std::endl;
    std::cout << std::endl;
  }
  else {
    std::cerr << "Something wrong solving system! info = " << info << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  #endif

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "coeffs_interior_ =" << std::endl;
  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    std::cout << std::setw(12) << coeffs_interior_[ii];
  }
  std::cout << std::endl << std::endl;
  #endif

  delete [] pp;
  pp = nullptr;

  return true;
}

bool mtk::Grad1D::ComputeRationalBasisNullSpace(void) {

  /// 1. Create generator vector for the first approximation.

  mtk::Real* gg{}; // Generator vector for the first Vandermonde matrix.

  try {
    gg = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(gg, mtk::kZero, sizeof(gg[0])*num_bndy_coeffs_);

  #ifdef MTK_PRECISION_DOUBLE
  gg[1] = 1.0/2.0;
  #else
  gg[1] = 1.0f/2.0f;
  #endif
  for (auto ii = 2; ii < num_bndy_coeffs_; ++ii) {
    gg[ii] = gg[ii - 1] + mtk::kOne;
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "gg =" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    std::cout << std::setw(12) << gg[ii];
  }
  std::cout << std::endl << std::endl;
  #endif

  /// 2. Create Vandermonde matrix.

  bool tran{true}; // Should I transpose the Vandermonde matrix.

  mtk::DenseMatrix aa_west_t(gg, num_bndy_coeffs_, order_accuracy_ + 1, tran);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "aa_west_t =" << std::endl;
  std::cout << aa_west_t << std::endl;
  #endif

  /// 3. QR-factorize the Vandermonde matrix.

  mtk::DenseMatrix qq_t(mtk::LAPACKAdapter::QRFactorDenseMatrix(aa_west_t));

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "qq_t = " << std::endl;
  std::cout << qq_t << std::endl;
  #endif

  /// 4.  Extract the basis for the null-space from Q matrix.

  int kk_num_rows{num_bndy_coeffs_};
  int kk_num_cols{dim_null_};

  mtk::DenseMatrix kk(kk_num_rows, kk_num_cols);

  // In the case of the gradient, even though we must solve for a null-space
  // of dimension 2, we must only extract ONE basis for the kernel.
  // We perform this extraction here:

  int aux_{kk_num_rows - kk_num_cols};
  for (auto ii = kk_num_rows - kk_num_cols; ii < kk_num_rows; ii++) {
    aux_--;
    for (auto jj = 0; jj < kk_num_rows; jj++) {
      kk.data()[jj*kk_num_cols + (kk_num_rows - kk_num_cols - aux_ - 1)] =
        qq_t.data()[ii*num_bndy_coeffs_ + jj];
    }
  }

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "kk =" << std::endl;
  std::cout << kk << std::endl;
  std::cout << "kk.num_rows() = " << kk.num_rows() << std::endl;
  std::cout << "kk.num_cols() = " << kk.num_cols() << std::endl;
  std::cout << std::endl;
  #endif

  /// 5. Scale null-space to make it rational.

  // Scale thus requesting that the last entries of the attained basis for the
  // null-space, adopt the pattern we require.
  // Essentially we will implement the following MATLAB pseudo-code:
  //  scalers = kk(num_bndy_approxs - (dim_null - 1):num_bndy_approxs,:)\B
  //  SK = kk*scalers
  // where SK is the scaled null-space.

  // In this point, we almost have all the data we need correctly allocated
  // in memory. We will create the matrix iden_, and elements we wish to scale
  // in the kk array. Using the concept of the leading dimension, we could just
  // use kk, with the correct leading dimension and that is it. BUT I DO NOT
  // GET how does it work. So I will just create a matrix with the content of
  // this array that we need, solve for the scalers and then scale the
  // whole kk:

  // We will then create memory for that sub-matrix of kk (subk).

  mtk::DenseMatrix subk(dim_null_, dim_null_);

  auto zz = 0;
  for (auto ii = order_accuracy_ + 1; ii < num_bndy_coeffs_; ii++) {
    for (auto jj = 0; jj < dim_null_; jj++) {
      subk.data()[zz*(dim_null_) + jj] = kk.data()[ii*(dim_null_) + jj];
    }
    zz++;
  }

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "subk =" << std::endl;
  std::cout << subk << std::endl;
  #endif

  subk.Transpose();

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "subk_t =" << std::endl;
  std::cout << subk << std::endl;
  #endif

  bool padded{false};
  tran = false;

  mtk::DenseMatrix iden(dim_null_, padded, tran);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "iden =" << std::endl;
  std::cout << iden << std::endl;
  #endif

  // Solve the system to compute the scalers.
  // An example of the system to solve, for k = 8, is:
  //
  // subk*scalers = iden or
  //
  // |  0.386018  -0.0339244   -0.129478 |           | 1 0 0 |
  // | -0.119774   0.0199423   0.0558632 |*scalers = | 0 1 0 |
  // | 0.0155708 -0.00349546 -0.00853182 |           | 0 0 1 |
  //
  // Notice this is a nrhs = 3 system.
  // Noteworthy: we do NOT ACTUALLY ALLOCATE space for the scalers... they
  // will be stored in the created identity matrix.
  // Let us first transpose subk (because of LAPACK):

  int info{mtk::LAPACKAdapter::SolveDenseSystem(subk, iden)};

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cout << "System successfully solved!" <<
      std::endl;
  } else {
    std::cerr << "Something went wrong solving system! info = " << info <<
      std::endl;
    std::cerr << "Exiting..." << std::endl;
    return false;
  }
  std::cout << std::endl;
  #endif

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Computed scalers:" << std::endl;
  std::cout << iden << std::endl;
  #endif

  // Multiply the two matrices to attain a scaled basis for null-space.

  rat_basis_null_space_ = mtk::BLASAdapter::RealDenseMM(kk, iden);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Rational basis for the null-space:" << std::endl;
  std::cout << rat_basis_null_space_ << std::endl;
  #endif

  // At this point, we have a rational basis for the null-space, with the
  // pattern we need! :)

  delete [] gg;
  gg = nullptr;

  return true;
}

bool mtk::Grad1D::ComputePreliminaryApproximations() {

  /// 1. Create generator vector for the first approximation.

  mtk::Real *gg{}; // Generator vector for the first approximation.

  try {
    gg = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(gg, mtk::kZero, sizeof(gg[0])*num_bndy_coeffs_);

  #ifdef MTK_PRECISION_DOUBLE
  gg[1] = 1.0/2.0;
  #else
  gg[1] = 1.0f/2.0f;
  #endif
  for (auto ii = 2; ii < num_bndy_coeffs_; ++ii) {
    gg[ii] = gg[ii - 1] + mtk::kOne;
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "gg0 =" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    std::cout << std::setw(12) << gg[ii];
  }
  std::cout << std::endl << std::endl;
  #endif

  // Allocate 2D array to store the collection of preliminary approximations.
  try {
    prem_apps_ = new mtk::Real[num_bndy_coeffs_*num_bndy_approxs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(prem_apps_,
         mtk::kZero,
         sizeof(prem_apps_[0])*num_bndy_coeffs_*num_bndy_approxs_);

  /// 2. Compute the dim_null near-the-boundary columns of the pi matrix.

  for (auto ll = 0; ll < num_bndy_approxs_; ++ll) {

    // Re-check new generator vector for every iteration except for the first.
    #if MTK_VERBOSE_LEVEL > 3
    if (ll > 0) {
      std::cout << "gg_" << ll << " =" << std::endl;
      for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
        std::cout << std::setw(12) << gg[ii];
      }
      std::cout << std::endl << std::endl;
    }
    #endif

    /// 3. Create the Vandermonde matrix for this iteration.

    bool transpose{false};

    mtk::DenseMatrix aa(gg,
                         num_bndy_coeffs_, order_accuracy_ + 1,
                         transpose);

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "aa_" << ll << " =" << std::endl;
    std::cout << aa << std::endl;
    #endif

    /// 4. New order-selector vector (gets re-written with LAPACK solutions).

    mtk::Real *ob{};

    auto ob_ld = num_bndy_coeffs_;

    try {
      ob = new mtk::Real[ob_ld];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(ob, mtk::kZero, sizeof(ob[0])*ob_ld);

    ob[1] = mtk::kOne;

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "ob = " << std::endl << std::endl;
    for (auto ii = 0; ii < ob_ld; ++ii) {
      std::cout << std::setw(12) << ob[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 5. Solving TT*rr = ob yields the columns rr of the kk matrix.

    // However, this is an under-determined system of equations. So we can not
    // use the same LAPACK routine (dgesv_). We will instead use dgels_, through
    // our LAPACKAdapter class.

    int info_{
      mtk::LAPACKAdapter::SolveRectangularDenseSystem(aa, ob, ob_ld)};

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!info_) {
      std::cout << "System successfully solved!" << std::endl << std::endl;
    } else {
      std::cerr << "Error solving system! info = " << info_ << std::endl;
      return false;
    }
    #endif

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "ob =" << std::endl;
    for (auto ii = 0; ii < ob_ld; ++ii) {
      std::cout << std::setw(12) << ob[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 6. Scale the kk matrix to make it a rational basis for null-space.

    // This implies a DAXPY operation. However, we must construct the arguments
    // for this operation.

    /// 7. Extract the last dim_null values of the pre-scaled ob.
    // Save them into the ob_bottom array:

    Real *ob_bottom{}; // Bottom part of the attained kernel used to scale it.

    try {
      ob_bottom = new mtk::Real[dim_null_];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(ob_bottom, mtk::kZero, sizeof(ob_bottom[0])*dim_null_);

    for (auto ii = 0; ii < dim_null_; ++ii) {
      ob_bottom[(dim_null_ - 1) - ii] = ob[num_bndy_coeffs_ - ii - 1];
    }

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "ob_bottom =" << std::endl;
    for (auto ii = 0; ii < dim_null_; ++ii) {
      std::cout << std::setw(12) << ob_bottom[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 8. Once we posses the bottom elements, we proceed with the scaling.

    // We must computed an scaled ob, sob, using the scaled null-space in
    // rat_basis_null_space_.
    // Such operation is: sob = ob - rat_basis_null_space_*ob_bottom
    // or:                 ob = -1.0*rat_basis_null_space_*ob_bottom + 1.0*ob
    // thus:                Y =    a*A    *x         +   b*Y (DAXPY).

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "Rational basis for the null-space:" << std::endl;
    std::cout << rat_basis_null_space_ << std::endl;
    #endif

    mtk::Real alpha{-mtk::kOne};
    mtk::Real beta{mtk::kOne};

    mtk::BLASAdapter::RealDenseMV(alpha, rat_basis_null_space_,
                                  ob_bottom, beta, ob);

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "scaled ob:" << std::endl;
    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      std::cout << std::setw(12) << ob[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    // We save the recently scaled solution, into an array containing these.
    // We can NOT start building the pi matrix, simply because I want that part
    // to be separated since its construction depends on the algorithm we want
    // to implement.

    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      prem_apps_[ii*num_bndy_approxs_ + ll] = ob[ii];
    }

    // After the first iteration, simply shift the entries of the last
    // generator vector used:
    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      gg[ii]--;
    }

    // Garbage collection for this loop:
    delete[] ob;
    ob = nullptr;

    delete[] ob_bottom;
    ob_bottom = nullptr;
  } // End of: for (ll = 0; ll < dim_null; ll++);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Matrix post-scaled preliminary apps: " << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < num_bndy_approxs_; ++jj) {
      std::cout << std::setw(12) << prem_apps_[ii*num_bndy_approxs_ + jj];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  #endif

  delete[] gg;
  gg = nullptr;

  return true;
}

bool mtk::Grad1D::ComputeWeights() {

  // Matrix to compute the weights as in the CRSA.
  mtk::DenseMatrix pi(num_bndy_coeffs_, num_bndy_coeffs_ - 1);

  /// 1. Construct the \f$ \mathbf{\Pi}\f$ matrix.

  // Assemble the pi matrix using:
  // 1. The collection of scaled preliminary approximations.
  // 2. The collection of coefficients approximating at the interior.
  // 3. The scaled basis for the null-space.

  // 1.1. Process array of scaled preliminary approximations.

  // These are queued in scaled_solutions. Each one of these, will be a column
  // of the pi matrix:
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < num_bndy_approxs_; ++jj) {
      pi.data()[ii*(2*(num_bndy_approxs_ - 1) + (order_accuracy_/2 + 1)) + jj] =
        prem_apps_[ii*num_bndy_approxs_ + jj];
    }
  }

  // 1.2. Add columns from known stencil approximating at the interior.

  // However, these must be padded by zeros, according to their position in the
  // final pi matrix:
  auto mm = 1;
  for (auto jj = num_bndy_approxs_; jj < order_accuracy_; ++jj) {
    for (auto ii = 0; ii < order_accuracy_; ++ii) {
      auto de = (ii + mm)*(2*(num_bndy_approxs_ - 1) +
        (order_accuracy_/2 + 1)) + jj;
      pi.data()[de] = coeffs_interior_[ii];
    }
    ++mm;
  }

  rat_basis_null_space_.OrderColMajor();

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Rational basis for the null-space (col. major):" << std::endl;
  std::cout << rat_basis_null_space_ << std::endl;
  #endif

  // 1.3. Add final set of columns: rational basis for null-space.

  for (auto jj = dim_null_ + (order_accuracy_/2 + 1);
       jj < num_bndy_coeffs_ - 1; ++jj) {
    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      auto og =
        (jj - (dim_null_ + (order_accuracy_/2 + 1)))*num_bndy_coeffs_ + ii;
      auto de = ii*(2*dim_null_ + (order_accuracy_/2 + 1)) + jj;
      pi.data()[de] = rat_basis_null_space_.data()[og];
    }
  }

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "coeffs_interior_ =" << std::endl;
  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    std::cout << std::setw(12) << coeffs_interior_[ii];
  }
  std::cout << std::endl << std::endl;
  #endif

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Constructed pi matrix for CRS Algorithm: " << std::endl;
  std::cout << pi << std::endl;
  #endif

  /// 2. Use interior stencil to build proper RHS vector \f$ \mathbf{h} \f$.

  // This imposes the mimetic condition.

  mtk::Real *hh{};  // Right-hand side to compute weights in the C{R,B}SA.

  try {
    hh = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(hh, mtk::kZero, sizeof(hh[0])*num_bndy_coeffs_);

  hh[0] = -mtk::kOne;
  for (auto ii = (order_accuracy_/2 + 2 - 1); ii < num_bndy_coeffs_; ++ii) {
    auto aux_xx = mtk::kZero;
    for (auto jj = 0; jj < ((ii - (order_accuracy_/2 - 1)) - 1); ++jj) {
      aux_xx += coeffs_interior_[jj];
    }
    hh[ii] = -mtk::kOne*aux_xx;
  }

  /// 3. Get weights (as **CRSA**):\f$ \mathbf{\Pi}\mathbf{q} = \mathbf{h} \f$.

  // That is, we construct a system, to solve for the weights.

  // Once again we face the challenge of solving with LAPACK. However, for the
  // CRSA, this matrix PI is over-determined, since it has more rows than
  // unknowns. However, according to the theory, the solution to this system is
  // unique. We will use dgels_.

  try {
    weights_cbs_ = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(weights_cbs_, mtk::kZero, sizeof(weights_cbs_[0])*num_bndy_coeffs_);

  int weights_ld{pi.num_cols() + 1};

  // Preserve hh.
  std::copy(hh, hh + weights_ld, weights_cbs_);

  pi.Transpose();

  int info{
    mtk::LAPACKAdapter::SolveRectangularDenseSystem(pi,
                                                    weights_cbs_, weights_ld)
  };

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cout << "System successfully solved!" << std::endl << std::endl;
  } else {
    std::cerr << "Error solving system! info = " << info << std::endl;
    return false;
  }
  #endif

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "hh =" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    std::cout << std::setw(11) << hh[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  // Preserve the original weights for research.

  try {
    weights_crs_ = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(weights_crs_, mtk::kZero, sizeof(weights_crs_[0])*num_bndy_coeffs_);

  std::copy(weights_cbs_, weights_cbs_ + (weights_ld - 1), weights_crs_);

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "weights_CRSA + lambda =" << std::endl;
  for (auto ii = 0; ii < weights_ld - 1; ++ii) {
    std::cout << std::setw(12) << weights_crs_[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 5. If required order is greater than critical order, start the **CBSA**.

  if (order_accuracy_ >= mtk::kCriticalOrderAccuracyGrad) {

    /// 6. Create \f$ \mathbf{\Phi} \f$ matrix from \f$ \mathbf{\Pi} \f$.

    mtk::DenseMatrix phi(order_accuracy_ + 1, order_accuracy_);

    // 6.1. Insert preliminary approximations to first set of columns.

    for (auto ii = 0; ii < order_accuracy_ + 1; ++ii) {
      for (auto jj = 0; jj < num_bndy_approxs_; ++jj) {
        phi.data()[ii*(order_accuracy_) + jj] =
          prem_apps_[ii*num_bndy_approxs_ + jj];
      }
    }

    // 6.2. Skip a column and negate preliminary approximations.

    for (auto jj = 0; jj < order_accuracy_ + 1; jj++) {
      for (auto ii = 1; ii < num_bndy_approxs_; ii++) {
        auto de = (ii+ order_accuracy_ - num_bndy_approxs_+ jj*order_accuracy_);
        auto og = (num_bndy_approxs_ - ii + (jj)*num_bndy_approxs_);
        phi.data()[de] = -prem_apps_[og];
      }
    }

    // 6.3. Flip negative columns up-down.

    for (auto ii = 0; ii < order_accuracy_/2; ii++) {
      for (auto jj = num_bndy_approxs_ + 1; jj < order_accuracy_; jj++) {
        auto aux = phi.data()[ii*order_accuracy_ + jj];
        phi.data()[ii*order_accuracy_ + jj] =
          phi.data()[(order_accuracy_ - ii)*order_accuracy_ + jj];
        phi.data()[(order_accuracy_ - ii)*order_accuracy_ + jj] = aux;
      }
    }

    // 6.4. Insert stencil.

    auto mm = 0;
    for (auto jj = num_bndy_approxs_; jj < num_bndy_approxs_ +  1; jj++) {
      for (auto ii = 0; ii < order_accuracy_ + 1; ii++) {
        if (ii == 0) {
          phi.data()[jj] = 0.0;
        } else {
          phi.data()[(ii + mm)*order_accuracy_ + jj] = coeffs_interior_[ii - 1];
        }
      }
      mm++;
    }

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "phi =" << std::endl;
    std::cout << phi << std::endl;
    #endif

    /// 7. Prepare constraint vector as in the CBSA: \f$ \mathbf{\Lambda}\f$.

    mtk::Real *lamed{};  // Used to build big lambda.

    try {
      lamed = new mtk::Real[num_bndy_approxs_ - 1];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(lamed, mtk::kZero, sizeof(lamed[0])*(num_bndy_approxs_ - 1));

    for (auto ii = 0; ii < num_bndy_approxs_ - 1; ++ii) {
      lamed[ii] = hh[ii + order_accuracy_ + 1] ;
    }

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "lamed =" << std::endl;
    for (auto ii = 0; ii < num_bndy_approxs_ - 1; ++ii) {
      std::cout << std::setw(12) << lamed[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      mtk::Real temp = mtk::kZero;
      for(auto jj = 0; jj < num_bndy_approxs_ - 1; ++jj) {
        temp = temp +
          lamed[jj]*rat_basis_null_space_.data()[jj*num_bndy_coeffs_ + ii];
      }
      hh[ii] = hh[ii] - temp;
    }

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "big_lambda =" << std::endl;
    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      std::cout << std::setw(12) << hh[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 8. Brute force search through all the rows of the \f$\Phi\f$ matrix.

    #ifdef MTK_VERBOSE_WEIGHTS
    int copy_result{1};
    #else
    int copy_result{};
    #endif

    int minrow_{std::numeric_limits<int>::infinity()};

    mtk::Real norm{mtk::BLASAdapter::RealNRM2(weights_cbs_,order_accuracy_)};
    mtk::Real minnorm{std::numeric_limits<mtk::Real>::infinity()};

    mtk::Real normerr_; // Norm of the error for the solution on each row.

    #ifdef MTK_VERBOSE_WEIGHTS
    std::ofstream table("grad_1d_" + std::to_string(order_accuracy_) +
      "_weights.tex");

    table << "\\begin{tabular}[c]{c";
    for (int ii = 1; ii <= order_accuracy_; ++ii) {
      table << 'c';
    }
    table << ":c}\n\\toprule\nRow & ";
    for (int ii = 1; ii <= order_accuracy_; ++ii) {
      table << "$ q_{" + std::to_string(ii) + "}$ &";
    }
    table << " Relative error \\\\\n\\midrule\n";
    #endif

    for(auto row_= 0; row_ < order_accuracy_ + 1; ++row_) {
      normerr_ = mtk::GLPKAdapter::SolveSimplexAndCompare(phi.data(),
                                                          order_accuracy_ + 1,
                                                          order_accuracy_,
                                                          order_accuracy_,
                                                          hh,
                                                          weights_cbs_,
                                                          row_,
                                                          mimetic_threshold_,
                                                          copy_result);
      mtk::Real aux{normerr_/norm};

      #if MTK_VERBOSE_LEVEL > 2
      std::cout << "Relative norm: " << aux << " " << std::endl;
      std::cout << std::endl;
      #endif

      num_feasible_sols_ = num_feasible_sols_ +
        (int) (normerr_ != std::numeric_limits<mtk::Real>::infinity());

      if (aux < minnorm) {
        minnorm = aux;
        minrow_= row_;
      }

      #ifdef MTK_VERBOSE_WEIGHTS
      table << std::to_string(row_ + 1) << " & ";
      if (normerr_ != std::numeric_limits<mtk::Real>::infinity()) {
        for (int ii = 1; ii <= order_accuracy_; ++ii) {
          table << std::to_string(weights_cbs_[ii - 1]) + " & ";
        }
        table << std::to_string(aux) << " \\\\" << std::endl;
      } else {
        table << "\\multicolumn{" << std::to_string(order_accuracy_) <<
          "}{c}{$\\emptyset$} & ";
        table << " - \\\\" << std::endl;
      }
      #endif
    }

    #ifdef MTK_VERBOSE_WEIGHTS
    table << "\\midrule" << std::endl;
    table << "CRS weights:";
    for (int ii = 1; ii <= order_accuracy_; ++ii) {
      table << " & " << std::to_string(weights_crs_[ii - 1]);
    }
    table << " & - \\\\\n\\bottomrule\n\\end{tabular}" << std::endl;
    table.close();
    #endif

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "weights_CBSA + lambda (after brute force search):" <<
      std::endl;
    for (auto ii = 0; ii < num_bndy_coeffs_ - 1; ++ii) {
      std::cout << std::setw(12) << weights_cbs_[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 9. Apply solution found from brute force search.

    // After we know which row yields the smallest relative norm that row is
    // chosen to be the objective function and the result of the optimizer is
    // chosen to be the new weights_.

    #if MTK_VERBOSE_LEVEL > 2
    std::cout << "Minimum Relative Norm " << minnorm << " found at row " <<
      minrow_ + 1 << std::endl;
    std::cout << std::endl;
    #endif

    copy_result = 1;
    normerr_ = mtk::GLPKAdapter::SolveSimplexAndCompare(phi.data(),
                                                        order_accuracy_ + 1,
                                                        order_accuracy_,
                                                        order_accuracy_,
                                                        hh,
                                                        weights_cbs_,
                                                        minrow_,
                                                        mimetic_threshold_,
                                                        copy_result);
    mtk::Real aux_{normerr_/norm};
    #if MTK_VERBOSE_LEVEL > 2
    std::cout << "Relative norm: " << aux_ << std::endl;
    std::cout << std::endl;
    #endif

    delete [] lamed;
    lamed = nullptr;
  }

  delete [] hh;
  hh = nullptr;

  return true;
}

bool mtk::Grad1D::ComputeStencilBoundaryGrid(void) {

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "weights_* + lambda =" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_ - 1; ++ii) {
    std::cout << std::setw(12) << weights_cbs_[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 1. Collect lambda values.

  mtk::Real *lambda{}; // Collection of bottom values from weights_.

  try {
    lambda = new mtk::Real[dim_null_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(lambda, mtk::kZero, sizeof(lambda[0])*dim_null_);

  for (auto ii = 0; ii < dim_null_; ++ii) {
    lambda[ii] = weights_cbs_[order_accuracy_ + ii];
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "lambda =" << std::endl;
  for (auto ii = 0; ii < dim_null_; ++ii) {
    std::cout << std::setw(12) << lambda[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 2. Compute alpha values.

  mtk::Real *alpha{}; // Collection of alpha values.

  try {
    alpha = new mtk::Real[dim_null_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(alpha, mtk::kZero, sizeof(alpha[0])*dim_null_);

  for (auto ii = 0; ii < dim_null_; ++ii) {
    alpha[ii] = lambda[ii]/weights_cbs_[ii] ;
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "alpha =" << std::endl;
  for (auto ii = 0; ii < dim_null_; ++ii) {
    std::cout << std::setw(12) << alpha[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 3. Compute the mimetic boundary approximations.

  try {
    mim_bndy_ = new mtk::Real[num_bndy_coeffs_*num_bndy_approxs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(mim_bndy_,
         mtk::kZero,
         sizeof(mim_bndy_[0])*num_bndy_coeffs_*num_bndy_approxs_);

  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < (num_bndy_approxs_ - 1); ++jj) {
      mim_bndy_[ii*num_bndy_approxs_ + jj] =
        prem_apps_[ii*num_bndy_approxs_ + jj] +
        alpha[jj]*rat_basis_null_space_.data()[jj*num_bndy_coeffs_ + ii];
    }
  }

  for(auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    mim_bndy_[ii*num_bndy_approxs_ + (num_bndy_approxs_ - 1)] =
      prem_apps_[ii*num_bndy_approxs_ + (num_bndy_approxs_ - 1)];
  }

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "Collection of mimetic approximations:" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < num_bndy_approxs_; ++jj) {
      std::cout << std::setw(13) << mim_bndy_[ii*num_bndy_approxs_ + jj];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 4. Compute the row-wise sum to double-check that the operator is mimetic.

  for (auto ii = 0; ii < num_bndy_approxs_; ++ii) {
    sums_rows_mim_bndy_.push_back(mtk::kZero);




    for (auto jj = 0; jj < num_bndy_coeffs_; ++jj) {
      sums_rows_mim_bndy_[ii] += mim_bndy_[jj*num_bndy_approxs_ + ii];
    }
  }

    mimetic_measure_ = *std::max_element(sums_rows_mim_bndy_.begin(),
                                      sums_rows_mim_bndy_.end());

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "Row-wise sum of mimetic approximations:" << std::endl;
  for (auto ii = 0; ii < num_bndy_approxs_; ++ii) {
    std::cout << std::setw(13) << sums_rows_mim_bndy_[ii];
  }
  std::cout << std::endl;
  std::cout << std::endl;
  #endif

  delete[] lambda;
  lambda = nullptr;

  delete[] alpha;
  alpha = nullptr;

  return true;
}

bool mtk::Grad1D::AssembleOperator(void) {

  // The output array will have this form:
  // 1. The first entry of the array will contain the used order kk.
  // 2. The second entry of the array will contain the collection of
  // approximating coefficients for the interior of the grid.
  // 3. The third entry will contain a collection of weights.
  // 4. The next dim_null - 1 entries will contain the collections of
  // approximating coefficients for the west boundary of the grid.

  gradient_length_ = 1 + order_accuracy_ + order_accuracy_ +
    num_bndy_approxs_*num_bndy_coeffs_;

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "gradient_length_ = " << gradient_length_ << std::endl;
  #endif

  try {
    gradient_ = new mtk::Real[gradient_length_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(gradient_, mtk::kZero, sizeof(gradient_[0])*gradient_length_);

  /// 1. The first entry of the array will contain the order of accuracy.

  gradient_[0] = order_accuracy_;

  /// 2. The second entry of the array will contain the collection of
  /// approximating coefficients for the interior of the grid.

  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    gradient_[ii + 1] = coeffs_interior_[ii];
  }

  /// 3. The third entry will contain the collection of weights.

  for (auto ii = 0; ii < order_accuracy_; ++ii) {
    gradient_[(order_accuracy_ + 1) + ii] = weights_cbs_[ii];
  }

  /// 4. The next dim_null + 1 entries will contain the collections of
  /// approximating coefficients for the west boundary of the grid.

  int offset{2*order_accuracy_ + 1};

  int aux {}; // Auxiliary variable.

  if (order_accuracy_ > mtk::kDefaultOrderAccuracy) {
    for (auto ii = 0; ii < num_bndy_approxs_ ; ii++) {
      for (auto jj = 0; jj < num_bndy_coeffs_; jj++) {
        gradient_[offset + aux] = mim_bndy_[jj*num_bndy_approxs_ + ii];
        aux++;
      }
    }
  } else {
    gradient_[offset + 0] = prem_apps_[0];
    gradient_[offset + 1] = prem_apps_[1];
    gradient_[offset + 2] = prem_apps_[2];
  }

  #if MTK_VERBOSE_LEVEL > 1
  std::cout << "1D " << order_accuracy_ << "-order grad built!" << std::endl;
  std::cout << std::endl;
  #endif

  return true;
}
