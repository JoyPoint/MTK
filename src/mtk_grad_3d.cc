/*!
\file mtk_grad_3d.cc

\brief Implements the class Grad3D.

This class implements a 3D gradient operator, constructed using the
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
#include "mtk_grad_3d.h"

mtk::Grad3D::Grad3D():
  order_accuracy_(),
  mimetic_threshold_() {}

mtk::Grad3D::Grad3D(const Grad3D &grad):
  order_accuracy_(grad.order_accuracy_),
  mimetic_threshold_(grad.mimetic_threshold_) {}

mtk::Grad3D::~Grad3D() {}

bool mtk::Grad3D::ConstructGrad3D(const mtk::UniStgGrid3D &grid,
                                  int order_accuracy,
                                  mtk::Real mimetic_threshold) {

  int num_cells_x = grid.num_cells_x();
  int num_cells_y = grid.num_cells_y();
  int num_cells_z = grid.num_cells_z();

  int mx = num_cells_x + 1;  // Gx vertical dimension.
  int nx = num_cells_x + 2;  // Gx horizontal dimension.
  int my = num_cells_y + 1;  // Gy vertical dimension.
  int ny = num_cells_y + 2;  // Gy horizontal dimension.
  int mz = num_cells_z + 1;  // Gz vertical dimension.
  int nz = num_cells_z + 2;  // Gz horizontal dimension.

  mtk::Grad1D grad;

  bool info = grad.ConstructGrad1D(order_accuracy, mimetic_threshold);

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cerr << "Mimetic grad could not be built." << std::endl;
    return info;
  }
  #endif

  auto west = grid.west_bndy();
  auto east = grid.east_bndy();
  auto south = grid.south_bndy();
  auto north = grid.east_bndy();
  auto bottom = grid.bottom_bndy();
  auto top = grid.top_bndy();

  mtk::UniStgGrid1D grid_x(west, east, num_cells_x);
  mtk::UniStgGrid1D grid_y(south, north, num_cells_y);
  mtk::UniStgGrid1D grid_z(bottom, top, num_cells_z);

  mtk::DenseMatrix Gx(grad.ReturnAsDenseMatrix(grid_x));
  mtk::DenseMatrix Gy(grad.ReturnAsDenseMatrix(grid_y));
  mtk::DenseMatrix Gz(grad.ReturnAsDenseMatrix(grid_z));

  bool padded{true};
  bool transpose{true};

  mtk::DenseMatrix tix(num_cells_x, padded, transpose);
  mtk::DenseMatrix tiy(num_cells_y, padded, transpose);
  mtk::DenseMatrix tiz(num_cells_z, padded, transpose);

  /// 1. Build preliminary staggering through the x direction.

  mtk::DenseMatrix aux1(mtk::DenseMatrix::Kron(tiz, tiy));
  mtk::DenseMatrix gx(mtk::DenseMatrix::Kron(aux1, Gx));

  /// 2. Build preliminary staggering through the y direction.

  mtk::DenseMatrix aux2(mtk::DenseMatrix::Kron(tiz, Gy));
  mtk::DenseMatrix gy(mtk::DenseMatrix::Kron(aux2, tix));

  /// 3. Build preliminary staggering through the z direction.

  mtk::DenseMatrix aux3(mtk::DenseMatrix::Kron(Gz, tiy));
  mtk::DenseMatrix gz(mtk::DenseMatrix::Kron(aux3, tix));

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "Gx: " << mx << " by " << nx << std::endl;
  std::cout << "Transpose Ix: " << num_cells_x << " by " << nx  << std::endl;
  std::cout << "Gy: " << my << " by " << ny << std::endl;
  std::cout << "Transpose Iy: " << num_cells_y << " by " << ny  << std::endl;
  std::cout << "Gz: " << mz << " by " << nz << std::endl;
  std::cout << "Transpose Iz: " << num_cells_z << " by " << nz  << std::endl;
  #endif

  /// 4. Actual operator: GG_xyz = [gx; gy; gz].

  int total_rows{mx*num_cells_y*num_cells_z +
                 num_cells_x*my*num_cells_z +
                 num_cells_x*num_cells_y*mz};
  int total_cols{nx*ny*nz};

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "Grad 3D: " << total_rows << " by " << total_cols << std::endl;
  #endif

  mtk::DenseMatrix g3d(total_rows, total_cols);

  for(auto ii = 0; ii < nx*ny*nz; ii++) {
    for(auto jj = 0; jj < mx*num_cells_y*num_cells_z; jj++) {
      g3d.SetValue(jj,ii, gx.GetValue(jj,ii));
    }

    int offset = mx*num_cells_y*num_cells_z;

    for(auto kk = 0; kk < num_cells_x*my*num_cells_z; kk++) {
      g3d.SetValue(kk + offset, ii, gy.GetValue(kk,ii));
    }

    offset += num_cells_x*my*num_cells_z;

    for(auto ll = 0; ll < num_cells_x*num_cells_y*mz; ll++) {
      g3d.SetValue(ll + offset, ii, gz.GetValue(ll,ii));
    }
  }

  gradient_ = g3d;

  return info;
}

mtk::DenseMatrix mtk::Grad3D::ReturnAsDenseMatrix() const {

  return gradient_;
}
