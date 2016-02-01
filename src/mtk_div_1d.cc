/*!
\file mtk_div_1d.cc

\brief Implements the class Div1D.

This class implements a 1D divergence matrix operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\todo Overload ostream operator as in mtk::Lap1D.

\todo Implement creation of \f$ \mathbf{\Lambda}\f$ w. mtk::BLASAdapter.
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

#include "mtk_div_1d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::Div1D &in) {

  int output_precision{4};
  int output_width{8};

  /// 1. Print order of accuracy.

  stream << "Order of accuracy: " << in.divergence_[0] << std::endl;

  /// 2. Print approximating coefficients for the interior.

  stream << "Interior stencil: " << std::endl;
  for (auto ii = 1; ii <= in.order_accuracy_; ++ii) {
    stream << std::setprecision(output_precision) << std::setw(output_width) <<
      in.divergence_[ii] << ' ';
  }
  stream << std::endl;

  if (in.order_accuracy_ > 2) {

    /// 3. Print mimetic weights.

    stream << "Weights:" << std::endl;
    for (auto ii = in.order_accuracy_ + 1; ii <= 2*in.order_accuracy_; ++ii) {
      stream << std::setprecision(output_precision) <<
        std::setw(output_width) << in.divergence_[ii] << ' ';
    }
    stream << std::endl;

    /// 4. Print mimetic approximations at the boundary.

    auto offset = (2*in.order_accuracy_ + 1);
    int mm{};
    for (auto ii = 0; ii < in.dim_null_; ++ii) {
      stream << "Mimetic boundary row " << ii + 1 << ":" << std::endl;
      for (auto jj = 0; jj < in.num_bndy_coeffs_; ++jj) {
        auto value = in.divergence_[offset + mm];
        stream << std::setprecision(output_precision) <<
          std::setw(output_width) << value << ' ';
        ++mm;
      }
      stream << std::endl;
      stream << "Sum of elements in row " << ii + 1 << ": " <<
        in.sums_rows_mim_bndy_[ii];
      stream << std::endl;
    }
  }

  return stream;
}
}

mtk::Div1D::Div1D():
  order_accuracy_(mtk::kDefaultOrderAccuracy),
  dim_null_(),
  num_bndy_coeffs_(),
  divergence_length_(),
  minrow_(),
  row_(),
  coeffs_interior_(),
  prem_apps_(),
  weights_crs_(),
  weights_cbs_(),
  mim_bndy_(),
  divergence_(),
  sums_rows_mim_bndy_(),
  mimetic_threshold_(mtk::kDefaultMimeticThreshold) {}

mtk::Div1D::Div1D(const Div1D &div):
  order_accuracy_(div.order_accuracy_),
  dim_null_(div.dim_null_),
  num_bndy_coeffs_(div.num_bndy_coeffs_),
  divergence_length_(div.divergence_length_),
  minrow_(div.minrow_),
  row_(div.row_),
  coeffs_interior_(div.coeffs_interior_),
  prem_apps_(div.prem_apps_),
  weights_crs_(div.weights_crs_),
  weights_cbs_(div.weights_cbs_),
  mim_bndy_(div.mim_bndy_),
  divergence_(div.divergence_),
  sums_rows_mim_bndy_(div.sums_rows_mim_bndy_),
  mimetic_threshold_(div.mimetic_threshold_) {}

mtk::Div1D::~Div1D() {

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

  delete[] divergence_;
  divergence_ = nullptr;
}

bool mtk::Div1D::ConstructDiv1D(int order_accuracy,
                                mtk::Real mimetic_threshold) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(order_accuracy < 2, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent((order_accuracy%2) != 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(mimetic_threshold <= mtk::kZero,
                      __FILE__, __LINE__, __func__);

  if (order_accuracy >= mtk::kCriticalOrderAccuracyDiv) {
    std::cout << "WARNING: Numerical accuracy is critical." << std::endl;
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

  // It is noteworthy, that the 2nd-order-accurate divergence operator has NO
  // approximation at the boundary, thus it has no weights. For this case, the
  // dimension of the null-space of the Vandermonde matrices used to compute the
  // approximating coefficients at the boundary is 0. Ergo, we compute this
  // number first and then decide if we must compute anything at the boundary.

  dim_null_ = order_accuracy_/2 - 1;

  if (dim_null_ > 0) {

    #ifdef MTK_PRECISION_DOUBLE
    num_bndy_coeffs_ = (int) (3.0*((mtk::Real) order_accuracy_)/2.0);
    #else
    num_bndy_coeffs_ = (int) (3.0f*((mtk::Real) order_accuracy_)/2.0f);
    #endif

    /// 2. Compute a rational basis for the null-space for the first matrix.

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

    abort_construction = ComputeRationalBasisNullSpace();

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!abort_construction) {
      std::cerr << "Could NOT complete stage 2.1." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      return false;
    }
    #endif

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

    abort_construction = ComputeStencilBoundaryGrid();

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!abort_construction) {
      std::cerr << "Could NOT complete stage 2.4." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      return false;
    }
    #endif

  } // End of: if (dim_null_ > 0);

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

int mtk::Div1D::num_bndy_coeffs() const {

  return num_bndy_coeffs_;
}

mtk::Real *mtk::Div1D::coeffs_interior() const {

  return coeffs_interior_;
}

mtk::Real *mtk::Div1D::weights_crs() const {

  return weights_crs_;
}

mtk::Real *mtk::Div1D::weights_cbs() const {

  return weights_cbs_;
}

mtk::DenseMatrix mtk::Div1D::mim_bndy() const {

  mtk::DenseMatrix xx(dim_null_, 3*order_accuracy_/2);

  auto counter = 0;
  for (auto ii = 0; ii < dim_null_; ++ii) {
    for(auto jj = 0; jj < 3*order_accuracy_/2; ++jj) {
      xx.SetValue(ii,jj, divergence_[2*order_accuracy_ + 1 + counter]);
      counter++;
    }
  }

  return xx;
}

std::vector<mtk::Real> mtk::Div1D::sums_rows_mim_bndy() const {

  return sums_rows_mim_bndy_;
}

mtk::DenseMatrix mtk::Div1D::ReturnAsDenseMatrix(
  const UniStgGrid1D &grid) const {

  int nn{grid.num_cells_x()}; // Number of cells on the grid.

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 3*order_accuracy_ - 1, __FILE__, __LINE__, __func__);
  #endif

  mtk::Real inv_delta_x{mtk::kOne/grid.delta_x()};

  int dd_num_rows = nn + 2;
  int dd_num_cols = nn + 1;
  int elements_per_row = num_bndy_coeffs_;
  int num_extra_rows = dim_null_;

  // Output matrix featuring sizes for divergence operators.
  mtk::DenseMatrix out(dd_num_rows, dd_num_cols);

  /// 1. Insert mimetic boundary at the west.

  auto ee_index = 0;
  for (auto ii = 1; ii < num_extra_rows + 1; ii++) {
    auto cc = 0;
    for(auto jj = 0 ; jj < dd_num_rows; jj++) {
      if( cc >= elements_per_row) {
        out.SetValue(ii, jj, mtk::kZero);
      } else {
        out.SetValue(ii,jj, mim_bndy_[ee_index++]*inv_delta_x);
        cc++;
      }
    }
  }

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = num_extra_rows + 1;
       ii < dd_num_rows - num_extra_rows - 1; ii++) {
    auto jj = ii - num_extra_rows - 1;
    for (auto cc = 0; cc < order_accuracy_; cc++, jj++) {
      out.SetValue(ii, jj, coeffs_interior_[cc]*inv_delta_x);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the mimetic boundaries.

  ee_index = 0;
  for (auto ii = dd_num_rows - 2; ii >= dd_num_rows - num_extra_rows - 1; ii--)
{
    auto cc = 0;
    for (auto jj = dd_num_cols - 1; jj >= 0; jj--) {
      if( cc >= elements_per_row) {
        out.SetValue(ii,jj,0.0);
      } else {
        out.SetValue(ii,jj,-mim_bndy_[ee_index++]*inv_delta_x);
        cc++;
      }
     }
  }

  return out;
}

mtk::DenseMatrix mtk::Div1D::ReturnAsDimensionlessDenseMatrix(
  int num_cells_x) const {

  int nn{num_cells_x}; // Number of cells on the grid.

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(nn < 3*order_accuracy_ - 1, __FILE__, __LINE__, __func__);
  #endif

  int dd_num_rows = nn + 2;
  int dd_num_cols = nn + 1;
  int elements_per_row = num_bndy_coeffs_;
  int num_extra_rows = dim_null_;

  // Output matrix featuring sizes for gradient operators.
  mtk::DenseMatrix out(dd_num_rows, dd_num_cols);

  /// 1. Insert mimetic boundary at the west.

  auto ee_index = 0;
  for (auto ii = 1; ii < num_extra_rows + 1; ii++) {
    auto cc = 0;
    for(auto jj = 0 ; jj < dd_num_rows; jj++) {
      if( cc >= elements_per_row) {
        out.SetValue(ii, jj, mtk::kZero);
      } else {
        out.SetValue(ii,jj, mim_bndy_[ee_index++]);
        cc++;
      }
    }
  }

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = num_extra_rows + 1;
       ii < dd_num_rows - num_extra_rows - 1; ii++) {
    auto jj = ii - num_extra_rows - 1;
    for (auto cc = 0; cc < order_accuracy_; cc++, jj++) {
      out.SetValue(ii, jj, coeffs_interior_[cc]);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the mimetic boundaries.

  ee_index = 0;
  for (auto ii = dd_num_rows - 2; ii >= dd_num_rows - num_extra_rows - 1; ii--)
  {
    auto cc = 0;
    for (auto jj = dd_num_cols - 1; jj >= 0; jj--) {
      if( cc >= elements_per_row) {
        out.SetValue(ii,jj,0.0);
      } else {
        out.SetValue(ii,jj,-mim_bndy_[ee_index++]);
        cc++;
      }
     }
  }

  return out;
}

bool mtk::Div1D::ComputeStencilInteriorGrid() {

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

  mtk::DenseMatrix vander_matrix(pp,
                                 order_accuracy_,
                                 order_accuracy_,
                                 transpose);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "vander_matrix = " << std::endl;
  std::cout << vander_matrix << std::endl;
  #endif

  /// 3. Create order-selector vector.

  try {
    coeffs_interior_ = new mtk::Real[order_accuracy_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(coeffs_interior_,
         mtk::kZero,
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

bool mtk::Div1D::ComputeRationalBasisNullSpace(void) {

  mtk::Real* gg{}; // Generator vector for the first Vandermonde matrix.

  /// 1. Create generator vector for the first approximation.

  try {
    gg = new mtk::Real[num_bndy_coeffs_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(gg, mtk::kZero, sizeof(gg[0])*num_bndy_coeffs_);

  #ifdef MTK_PRECISION_DOUBLE
  gg[0] = -1.0/2.0;
  #else
  gg[0] = -1.0f/2.0f;
  #endif
  for (auto ii = 1; ii < num_bndy_coeffs_; ++ii) {
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

  mtk::DenseMatrix vv_west_t(gg, num_bndy_coeffs_, order_accuracy_ + 1, tran);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "vv_west_t =" << std::endl;
  std::cout << vv_west_t << std::endl;
  #endif

  /// 3. QR-factorize the Vandermonde matrix.

  mtk::DenseMatrix qq_t(mtk::LAPACKAdapter::QRFactorDenseMatrix(vv_west_t));

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "QQ^T = " << std::endl;
  std::cout << qq_t << std::endl;
  #endif

  /// 4.  Extract the basis for the null-space from Q matrix.

  int KK_num_rows_{num_bndy_coeffs_};
  int KK_num_cols_{dim_null_};

  mtk::DenseMatrix KK(KK_num_rows_, KK_num_cols_);

  for (auto ii = num_bndy_coeffs_ - dim_null_; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < num_bndy_coeffs_; ++jj) {
      KK.data()[jj*dim_null_ + (ii - (num_bndy_coeffs_ - dim_null_))] =
          qq_t.data()[ii*num_bndy_coeffs_ + jj];
    }
  }

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "KK =" << std::endl;
  std::cout << KK << std::endl;
  std::cout << "KK.num_rows() = " << KK.num_rows() << std::endl;
  std::cout << "KK.num_cols() = " << KK.num_cols() << std::endl;
  std::cout << std::endl;
  #endif

  /// 5. Scale null-space to make it rational.

  // Scale thus requesting that the last entries of the attained basis for the
  // null-space, adopt the pattern we require.
  // Essentially we will implement the following MATLAB pseudo-code:
  //  scalers = KK(num_bndy_approxs - (dim_null - 1):num_bndy_approxs,:)\B
  //  SK = KK*scalers
  // where SK is the scaled null-space.

  // In this point, we almost have all the data we need correctly allocated
  // in memory. We will create the matrix II_, and elements we wish to scale in
  // the KK array. Using the concept of the leading dimension, we could just
  // use KK, with the correct leading dimension and that is it. BUT I DO NOT
  // GET how does it work. So I will just create a matrix with the content of
  // this array that we need, solve for the scalers and then scale the
  // whole KK:

  // We will then create memory for that sub-matrix of KK (SUBK).

  mtk::DenseMatrix SUBK(dim_null_,dim_null_);

  for (auto ii = num_bndy_coeffs_ - dim_null_; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < dim_null_; ++jj) {
      SUBK.data()[(ii - (num_bndy_coeffs_ - dim_null_))*dim_null_ + jj] =
          KK.data()[ii*dim_null_ + jj];
    }
  }

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "SUBK =" << std::endl;
  std::cout << SUBK << std::endl;
  #endif

  SUBK.Transpose();

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "SUBK^T =" << std::endl;
  std::cout << SUBK << std::endl;
  #endif

  bool padded{false};
  tran = false;

  mtk::DenseMatrix II(dim_null_, padded, tran);

  #if MTK_VERBOSE_LEVEL > 4
  std::cout << "II =" << std::endl;
  std::cout << II << std::endl;
  #endif

  // Solve the system to compute the scalers.
  // An example of the system to solve, for k = 8, is:
  //
  // SUBK*scalers = II_ or
  //
  // |  0.386018  -0.0339244   -0.129478 |           | 1 0 0 |
  // | -0.119774   0.0199423   0.0558632 |*scalers = | 0 1 0 |
  // | 0.0155708 -0.00349546 -0.00853182 |           | 0 0 1 |
  //
  // Notice this is a nrhs = 3 system.
  // Noteworthy: we do NOT ACTUALLY ALLOCATE space for the scalers... they
  // will be stored in the created identity matrix.
  // Let us first transpose SUBK (because of LAPACK):

  int info{mtk::LAPACKAdapter::SolveDenseSystem(SUBK, II)};

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
  std::cout << II << std::endl;
  #endif

  // Multiply the two matrices to attain a scaled basis for null-space.

  rat_basis_null_space_ = mtk::BLASAdapter::RealDenseMM(KK, II);

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

bool mtk::Div1D::ComputePreliminaryApproximations(void) {

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
  gg[0] = -1.0/2.0;
  #else
  gg[0] = -1.0f/2.0f;
  #endif
  for (auto ii = 1; ii < num_bndy_coeffs_; ++ii) {
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
    prem_apps_ = new mtk::Real[num_bndy_coeffs_*dim_null_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(prem_apps_,
         mtk::kZero,
         sizeof(prem_apps_[0])*num_bndy_coeffs_*dim_null_);

  /// 2. Compute the dim_null near-the-boundary columns of the pi matrix.

  for (auto ll = 0; ll < dim_null_; ++ll) {

    // Re-check new generator vector for every iteration except for the first.
    #if MTK_VERBOSE_LEVEL > 3
    if (ll > 0) {
      std::cout << "gg" << ll << " =" << std::endl;
      for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
        std::cout << std::setw(12) << gg[ii];
      }
      std::cout << std::endl << std::endl;
    }
    #endif

    /// 3. Create the Vandermonde matrix for this iteration.

    bool transpose{false};

    mtk::DenseMatrix AA_(gg,
                         num_bndy_coeffs_, order_accuracy_ + 1,
                         transpose);

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "AA_" << ll << " =" << std::endl;
    std::cout << AA_ << std::endl;
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

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "ob = " << std::endl << std::endl;
    for (auto ii = 0; ii < ob_ld; ++ii) {
      std::cout << std::setw(12) << ob[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 5. Solving TT*rr = ob yields the columns rr of the KK matrix.

    // However, this is an under-determined system of equations. So we can not
    // use the same LAPACK routine (dgesv_). We will instead use dgels_, through
    // our LAPACKAdapter class.

    int info_{
      mtk::LAPACKAdapter::SolveRectangularDenseSystem(AA_, ob, ob_ld)};

    #ifdef MTK_PERFORM_PREVENTIONS
    if (!info_) {
      std::cout << "System successfully solved!" << std::endl << std::endl;
    } else {
      std::cerr << "Error solving system! info = " << info_ << std::endl;
    }
    #endif

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "ob =" << std::endl;
    for (auto ii = 0; ii < ob_ld; ++ii) {
      std::cout << std::setw(12) << ob[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    /// 6. Scale the KK matrix to make it a rational basis for null-space.

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

    #if MTK_VERBOSE_LEVEL > 3
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
      prem_apps_[ii*dim_null_ + ll] = ob[ii];
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
    for (auto jj = 0; jj < dim_null_; ++jj) {
      std::cout << std::setw(12) << prem_apps_[ii*dim_null_ + jj];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  #endif

  delete[] gg;
  gg = nullptr;

  return true;
}

bool mtk::Div1D::ComputeWeights(void) {

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
    for (auto jj = 0; jj < dim_null_; ++jj) {
      pi.data()[ii*(2*dim_null_ + (order_accuracy_/2 + 1)) + jj] =
        prem_apps_[ii*dim_null_ + jj];
    }
  }

  // 1.2. Add columns from known stencil approximating at the interior.

  // However, these must be padded by zeros, according to their position in the
  // final pi matrix:
  auto mm = 0;
  for (auto jj = dim_null_; jj < order_accuracy_; ++jj) {
    for (auto ii = 0; ii < order_accuracy_; ++ii) {
      pi.data()[(ii + mm)*(2*dim_null_ + (order_accuracy_/2 + 1)) + jj] =
        coeffs_interior_[ii];
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
       jj < num_bndy_coeffs_ - 1;
       ++jj) {
    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      auto og =
        (jj - (dim_null_ + (order_accuracy_/2 + 1)))*num_bndy_coeffs_ + ii;
      auto de = ii*(2*dim_null_ + (order_accuracy_/2 + 1)) + jj;
      pi.data()[de] = rat_basis_null_space_.data()[og];
    }
  }

  #if MTK_VERBOSE_LEVEL > 3
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

  int info{mtk::LAPACKAdapter::SolveRectangularDenseSystem(pi,
                                                           weights_cbs_,
                                                           weights_ld)};

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cout << "System successfully solved!" << std::endl << std::endl;
  } else {
    std::cerr << "Error solving system! info = " << info << std::endl;
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

  if (order_accuracy_ >= mtk::kCriticalOrderAccuracyDiv) {

    /// 6. Create \f$ \mathbf{\Phi} \f$ matrix from \f$ \mathbf{\Pi} \f$.

    mtk::DenseMatrix phi(order_accuracy_ + 1, order_accuracy_);

    for (auto ii = 0; ii < order_accuracy_ + 1; ++ii) {
      for (auto jj = 0; jj < dim_null_; ++jj) {
        phi.data()[ii*(order_accuracy_) + jj] = prem_apps_[ii*dim_null_ + jj];
      }
    }

    int aux{};  // Auxiliary variable.
    for (auto jj = dim_null_; jj < dim_null_ + 2; ++jj) {
      for (auto ii = 0; ii < order_accuracy_; ++ii) {
        phi.data()[(ii + aux)*order_accuracy_ + jj] = coeffs_interior_[ii];
      }
      ++aux;
    }

    for(auto jj=order_accuracy_ - 1; jj >=order_accuracy_ - dim_null_; jj--) {
      for(auto ii=0; ii<order_accuracy_ + 1; ++ii) {
        phi.data()[ii*order_accuracy_+jj] = mtk::kZero;
      }
    }

    for (auto jj = 0; jj < order_accuracy_ + 1; ++jj) {
      for (auto ii = 0; ii < dim_null_; ++ii) {
        phi.data()[(ii + order_accuracy_ - dim_null_ + jj*order_accuracy_)] =
          -prem_apps_[(dim_null_ - ii - 1 + jj*dim_null_)];
      }
    }

    for(auto ii = 0; ii < order_accuracy_/2; ++ii) {
      for (auto jj = dim_null_ + 2; jj < order_accuracy_; ++jj) {
        auto swap = phi.data()[ii*order_accuracy_+jj];
        phi.data()[ii*order_accuracy_ + jj] =
          phi.data()[(order_accuracy_-ii)*order_accuracy_+jj];
        phi.data()[(order_accuracy_-ii)*order_accuracy_+jj] = swap;
      }
    }

    #if MTK_VERBOSE_LEVEL > 4
    std::cout << "Constructed PHI matrix for CBS Algorithm: " << std::endl;
    std::cout << phi << std::endl;
    #endif

    /// 7. Prepare constraint vector as in the CBSA: \f$ \mathbf{\Lambda}\f$.

    mtk::Real *lamed{};  // Used to build big lambda.

    try {
      lamed = new mtk::Real[dim_null_];
    } catch (std::bad_alloc &memory_allocation_exception) {
      std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
        std::endl;
      std::cerr << memory_allocation_exception.what() << std::endl;
    }
    memset(lamed, mtk::kZero, sizeof(lamed[0])*dim_null_);

    for (auto ii = 0; ii < dim_null_; ++ii) {
      lamed[ii] = hh[ii + order_accuracy_ + 1] ;
    }

    #if MTK_VERBOSE_LEVEL > 3
    std::cout << "lamed =" << std::endl;
    for (auto ii = 0; ii < dim_null_; ++ii) {
      std::cout << std::setw(12) << lamed[ii] << std::endl;
    }
    std::cout << std::endl;
    #endif

    for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
      mtk::Real temp = mtk::kZero;
      for(auto jj = 0; jj < dim_null_; ++jj) {
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

    #ifdef MTK_VERBOSE_WEIGHTS
    int copy_result{1};
    #else
    int copy_result{};
    #endif

    mtk::Real normerr_; // Norm of the error for the solution on each row.

    /// 8. Brute force search through all the rows of the \f$\Phi\f$ matrix.

    int minrow_{std::numeric_limits<int>::infinity()};

    mtk::Real norm_{mtk::BLASAdapter::RealNRM2(weights_crs_,order_accuracy_)};
    mtk::Real minnorm_{std::numeric_limits<mtk::Real>::infinity()};

    #ifdef MTK_VERBOSE_WEIGHTS
    std::ofstream table("div_1d_" + std::to_string(order_accuracy_) +
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
      mtk::Real aux{normerr_/norm_};

      #if MTK_VERBOSE_LEVEL > 2
      std::cout << "Relative norm: " << aux << " " << std::endl;
      std::cout << std::endl;
      #endif

      if (aux < minnorm_) {
        minnorm_ = aux;
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
    std::cout << "Minimum Relative Norm " << minnorm_ << " found at row " <<
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
    mtk::Real aux_{normerr_/norm_};
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

bool mtk::Div1D::ComputeStencilBoundaryGrid(void) {

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "weights_CBSA + lambda =" << std::endl;
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
    mim_bndy_ = new mtk::Real[num_bndy_coeffs_*dim_null_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(mim_bndy_,
         mtk::kZero,
         sizeof(mim_bndy_[0])*num_bndy_coeffs_*dim_null_);

  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < dim_null_; ++jj) {
      mim_bndy_[ii*dim_null_ + jj] =
        prem_apps_[ii*dim_null_ + jj] +
        alpha[jj]*rat_basis_null_space_.data()[jj*num_bndy_coeffs_ + ii];
    }
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "Collection of mimetic approximations:" << std::endl;
  for (auto ii = 0; ii < num_bndy_coeffs_; ++ii) {
    for (auto jj = 0; jj < dim_null_; ++jj) {
      std::cout << std::setw(13) << mim_bndy_[ii*dim_null_ + jj];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 4. Compute the row-wise sum to double-check the operator is mimetic.

  for (auto ii = 0; ii < dim_null_; ++ii) {
    sums_rows_mim_bndy_.push_back(mtk::kZero);
    for (auto jj = 0; jj < num_bndy_coeffs_; ++jj) {
      sums_rows_mim_bndy_[ii] += mim_bndy_[jj*dim_null_ + ii];
    }
  }

  #if MTK_VERBOSE_LEVEL > 3
  std::cout << "Row-wise sum of mimetic approximations:" << std::endl;
  for (auto ii = 0; ii < dim_null_; ++ii) {
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

bool mtk::Div1D::AssembleOperator(void) {

  // The output array will have this form:
  // 1. The first entry of the array will contain used order order_accuracy_.
  // 2. The second entry of the array will contain the collection of
  // approximating coefficients for the interior of the grid.
  // 3. IF order_accuracy_ > 2, then the third entry will contain a collection
  // of weights.
  // 4. IF order_accuracy_ > 2, the next dim_null_ entries will contain the
  // collections of approximating coefficients for the west boundary of the
  // grid.

  if (order_accuracy_ > mtk::kDefaultOrderAccuracy) {
    divergence_length_ =
      1 + order_accuracy_ + order_accuracy_ + dim_null_*num_bndy_coeffs_;
  } else {
    divergence_length_ = 1 + order_accuracy_;
  }

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "divergence_length_ = " << divergence_length_ << std::endl;
  std::cout << std::endl;
  #endif

  try {
    divergence_ = new double[divergence_length_];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(divergence_, mtk::kZero, sizeof(divergence_[0])*divergence_length_);

  /// 1. The first entry of the array will contain the order of accuracy.

  divergence_[0] = order_accuracy_;

  /// 2. The second entry the collection of coefficients for interior of grid.

  for (auto ii = 0; ii < order_accuracy_; ++ii) {
  divergence_[ii + 1] = coeffs_interior_[ii];
  }

  /// 3. If order_accuracy_ > 2, then third entry is the collection of weights.

  if (order_accuracy_ > 2) {
    for (auto ii = 0; ii < order_accuracy_; ++ii) {
      divergence_[(1 + order_accuracy_) + ii] = weights_cbs_[ii];
    }
  }

  /// 4. If order_accuracy_ > 2, next dim_null_ entries is approximating
  /// coefficients for the west boundary of the grid.

  if (order_accuracy_ > 2) {
    auto offset = (2*order_accuracy_ + 1);
    int mm{};
    for (auto ii = 0; ii < dim_null_; ++ii) {
      for (auto jj = 0; jj < num_bndy_coeffs_; ++jj) {
        divergence_[offset + (mm)] = mim_bndy_[jj*dim_null_ + ii];
        ++mm;
      }
    }
  }

  #if MTK_VERBOSE_LEVEL > 1
  std::cout << "1D " << order_accuracy_ << "-order div built!" << std::endl;
  std::cout << std::endl;
  #endif

  return true;
}
