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

bool mtk::Grad2D::ConstructGrad2D(const mtk::UniStgGrid2D &grid,
                                  int order_accuracy,
                                  mtk::Real mimetic_threshold) {

  int num_cells_x = grid.num_cells_x();
  int num_cells_y = grid.num_cells_y();

  int mx = num_cells_x + 1;  // Gx vertical dimension
  int nx = num_cells_x + 2;  // Gx horizontal dimension
  int my = num_cells_y + 1;  // Gy vertical dimension
  int ny = num_cells_y + 2;  // Gy horizontal dimension

  mtk::Grad1D grad;

  bool info = grad.ConstructGrad1D(order_accuracy, mimetic_threshold);

  if (!info) {
    std::cerr << "Mimetic grad could not be built." << std::endl;
    return info;
  }

  auto west = grid.west_bndy();
  auto east = grid.east_bndy();
  auto south = grid.south_bndy();
  auto north = grid.east_bndy();

  mtk::UniStgGrid1D grid_x(west, east, num_cells_x);
  mtk::UniStgGrid1D grid_y(south, north, num_cells_y);

  mtk::DenseMatrix Gx(grad.ReturnAsDenseMatrix(grid_x));
  mtk::DenseMatrix Gy(grad.ReturnAsDenseMatrix(grid_y));

  bool padded{true};
  bool transpose{true};

  mtk::DenseMatrix tix(num_cells_x, padded, transpose);
  mtk::DenseMatrix tiy(num_cells_y, padded, transpose);

  mtk::DenseMatrix gxy(mtk::DenseMatrix::Kron(tiy, Gx));
  mtk::DenseMatrix gyx(mtk::DenseMatrix::Kron(Gy, tix));

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Gx :" << mx << "by " << nx << std::endl;
  std::cout << "Transpose Iy : " << num_cells_y<< " by " << ny  << std::endl;
  std::cout << "Gy :" << my << "by " << ny << std::endl;
  std::cout << "Transpose Ix : " << num_cells_x<< " by " << nx  << std::endl;
  std::cout << "Kronecker dimensions Grad 2D" <<
    mx*num_cells_y + my*num_cells_x << " by " <<  nx*ny <<std::endl;
  #endif

  mtk::DenseMatrix g2d(mx*num_cells_y + my*num_cells_x, nx*ny);

  for(auto ii = 0; ii < nx*ny; ii++) {
    for(auto jj = 0; jj < mx*num_cells_y; jj++) {
      g2d.SetValue(jj,ii, gxy.GetValue(jj,ii));
    }
    for(auto kk = 0; kk < my*num_cells_x; kk++) {
      g2d.SetValue(kk + mx*num_cells_y, ii, gyx.GetValue(kk,ii));
    }
  }

  gradient_ = g2d;

  return info;
}

mtk::DenseMatrix mtk::Grad2D::ReturnAsDenseMatrix() const {

  return gradient_;
}
