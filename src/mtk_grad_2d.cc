/*!
\file mtk_grad_2d.cc

\brief Implements the class Grad2D.

This class implements a 2D gradient operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_grad_1d.h"
#include "mtk_grad_2d.h"

mtk::Grad2D::Grad2D():
  order_accuracy_(),
  mimetic_threshold_() {}

mtk::Grad2D::Grad2D(const Grad2D &grad):
  order_accuracy_(grad.order_accuracy_),
  mimetic_threshold_(grad.mimetic_threshold_) {}

mtk::Grad2D::~Grad2D() {}

mtk::DenseMatrix mtk::Grad2D::ConstructGrad2D(const mtk::UniStgGrid2D &grid,
                                              int order_accuracy,
                                              mtk::Real mimetic_threshold) {

  int NumCellsX = grid.num_cells_x();
  int NumCellsY = grid.num_cells_y();

  int mx = NumCellsX + 1;  // Gx vertical dimension
  int nx = NumCellsX + 2;  // Gx horizontal dimension
  int my = NumCellsY + 1;  // Gy vertical dimension
  int ny = NumCellsY + 2;  // Gy horizontal dimension

  mtk::Grad1D grad;

  bool info = grad.ConstructGrad1D(order_accuracy, mimetic_threshold);

  if (!info) {
    std::cerr << "Mimetic grad could not be built." << std::endl;
  }

  auto West = grid.west_bndy_x();
  auto East = grid.east_bndy_x();
  auto South = grid.south_bndy_y();
  auto North = grid.east_bndy_x();

  mtk::UniStgGrid1D grid_x(West, East, NumCellsX);
  mtk::UniStgGrid1D grid_y(South, North, NumCellsY);

  mtk::DenseMatrix Gx(grad.ReturnAsDenseMatrix(grid_x));
  mtk::DenseMatrix Gy(grad.ReturnAsDenseMatrix(grid_y));

  bool padded{true};
  bool transpose{true};

  mtk::DenseMatrix TIx(NumCellsX, padded, transpose);
  mtk::DenseMatrix TIy(NumCellsY, padded, transpose);

  mtk::DenseMatrix Gxy(mtk::DenseMatrix::Kron(TIy, Gx));
  mtk::DenseMatrix Gyx(mtk::DenseMatrix::Kron(Gy, TIx));

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Gx :" << mx << "by " << nx << std::endl;
  std::cout << "Transpose Iy : " << NumCellsY<< " by " << ny  << std::endl;
  std::cout << "Gy :" << my << "by " << ny << std::endl;
  std::cout << "Transpose Ix : " << NumCellsX<< " by " << nx  << std::endl;
  std::cout << "Kronecker dimensions Grad 2D" <<
  mx*NumCellsY + my*NumCellsX << " by " <<  nx*ny <<std::endl;
  #endif

  mtk::DenseMatrix G2D(mx*NumCellsY + my*NumCellsX, nx*ny);

  for(auto ii = 0; ii < nx*ny; ii++) {
    for(auto jj = 0; jj < mx*NumCellsY; jj++) {
      G2D.SetValue(jj,ii, Gxy.GetValue(jj,ii));
    }
    for(auto kk = 0; kk < my*NumCellsX; kk++) {
      G2D.SetValue(kk + mx*NumCellsY, ii, Gyx.GetValue(kk,ii));
    }
  }

  gradient_ = G2D;

  return gradient_;
}

mtk::DenseMatrix mtk::Grad2D::ReturnAsDenseMatrix() {

  return gradient_;
}
