/*!
\file mtk_interp_1d.cc

\brief Includes the implementation of the class Interp1D.

This class implements a 1D interpolation operator.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\author: Johnny Corbino - jcorbino at mail dot sdsu dot edu
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

#include <cstring>

#include "mtk_tools.h"

#include "mtk_interp_1d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::Interp1D &in) {

  /// 1. Print approximating coefficients for the interior.

  stream << "coeffs_interior_[1:" << in.order_accuracy_ << "] = ";
  for (auto ii = 0; ii < in.order_accuracy_; ++ii) {
    stream << std::setw(9) << in.coeffs_interior_[ii] << " ";
  }
  stream << std::endl;

  return stream;
}
}

mtk::Interp1D::Interp1D():
  dir_interp_(mtk::SCALAR_TO_VECTOR),
  order_accuracy_(mtk::kDefaultOrderAccuracy),
  coeffs_interior_(nullptr) {}

mtk::Interp1D::Interp1D(const Interp1D &interp):
  dir_interp_(interp.dir_interp_),
  order_accuracy_(interp.order_accuracy_),
  coeffs_interior_(interp.coeffs_interior_) {}

mtk::Interp1D::~Interp1D() {

  delete[] coeffs_interior_;
  coeffs_interior_ = nullptr;
}

bool mtk::Interp1D::ConstructInterp1D(int order_accuracy, mtk::DirInterp dir) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(order_accuracy < 2, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent((order_accuracy%2) != 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(dir < mtk::SCALAR_TO_VECTOR &&
                      dir > mtk::VECTOR_TO_SCALAR,
                      __FILE__, __LINE__, __func__);

  std::cout << "order_accuracy_ = " << order_accuracy << std::endl;
  #endif

  order_accuracy_ = order_accuracy;

  /// 1. Compute stencil for the interior cells.

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

  for (int ii = 0; ii < order_accuracy_; ++ii) {
    coeffs_interior_[ii] = mtk::kOne;
  }

  return true;
}

mtk::Real *mtk::Interp1D::coeffs_interior() const {

  return coeffs_interior_;
}

mtk::DenseMatrix mtk::Interp1D::ReturnAsDenseMatrix(const UniStgGrid1D &grid) {

  int nn{grid.num_cells_x()}; // Number of cells on the grid.

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  #endif

  int gg_num_rows{};  // Number of rows.
  int gg_num_cols{};  // Number of columns.

  if (dir_interp_ == mtk::SCALAR_TO_VECTOR) {
    gg_num_rows = nn + 1;
    gg_num_cols = nn + 2;
  } else {
    gg_num_rows = nn + 2;
    gg_num_cols = nn + 1;
  }

  // Output matrix featuring sizes for gradient operators.
  
  mtk::DenseMatrix out(gg_num_rows, gg_num_cols);

  /// 1. Preserve values at the boundary.

  out.SetValue(0, 0, mtk::kOne);

  /// 2. Insert coefficients for the interior of the grid.

  for (auto ii = 1; ii < gg_num_rows - 1; ++ii) {
    for(auto jj = ii ; jj < order_accuracy_ + ii; ++jj) {
      out.SetValue(ii, jj, mtk::kOne/order_accuracy_);
    }
  }

  /// 3. Impose center-skew symmetry by permuting the boundaries.

  out.SetValue(gg_num_rows - 1, gg_num_cols - 1, mtk::kOne);

  return out;
}
