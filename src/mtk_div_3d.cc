/*!
\file mtk_div_3d.cc

\brief Implements the class Div3D.

This class implements a 3D divergence operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).

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

#include <iostream>
#include <iomanip>

#include "mtk_foundations.h"
#include "mtk_div_1d.h"
#include "mtk_div_3d.h"

mtk::Div3D::Div3D():
  order_accuracy_(),
  mimetic_threshold_() {}

mtk::Div3D::Div3D(const Div3D &grad):
  order_accuracy_(grad.order_accuracy_),
  mimetic_threshold_(grad.mimetic_threshold_) {}

mtk::Div3D::~Div3D() {}

bool mtk::Div3D::ConstructDiv3D(const mtk::UniStgGrid3D &grid,
                                int order_accuracy,
                                mtk::Real mimetic_threshold) {

  int num_cells_x = grid.num_cells_x();
  int num_cells_y = grid.num_cells_y();
  int num_cells_z = grid.num_cells_z();

  int mx = num_cells_x + 1;  // Dx vertical dimension.
  int nx = num_cells_x + 2;  // Dx horizontal dimension.
  int my = num_cells_y + 1;  // Dy vertical dimension.
  int ny = num_cells_y + 2;  // Dy horizontal dimension.
  int mz = num_cells_z + 1;  // Dz vertical dimension.
  int nz = num_cells_z + 2;  // Dz horizontal dimension.

  mtk::Div1D div;

  bool info = div.ConstructDiv1D(order_accuracy, mimetic_threshold);

  #ifdef MTK_PERFORM_PREVENTIONS
  if (!info) {
    std::cerr << "Mimetic div could not be built." << std::endl;
    return info;
  }
  #endif

  auto west = grid.west_bndy();
  auto east = grid.east_bndy();
  auto south = grid.south_bndy();
  auto north = grid.east_bndy();
  auto bottom = grid.bottom_bndy();
  auto top = grid.top_bndy();

  mtk::UniStgGrid1D grid_x(west, east, num_cells_x, mtk::FieldNature::VECTOR);
  mtk::UniStgGrid1D grid_y(south, north, num_cells_y, mtk::FieldNature::VECTOR);
  mtk::UniStgGrid1D grid_z(bottom, top, num_cells_z, mtk::FieldNature::VECTOR);

  mtk::DenseMatrix Dx(div.ReturnAsDenseMatrix(grid_x));
  mtk::DenseMatrix Dy(div.ReturnAsDenseMatrix(grid_y));
  mtk::DenseMatrix Dz(div.ReturnAsDenseMatrix(grid_z));

  bool padded{true};
  bool transpose{false};

  mtk::DenseMatrix ix(num_cells_x, padded, transpose);
  mtk::DenseMatrix iy(num_cells_y, padded, transpose);
  mtk::DenseMatrix iz(num_cells_z, padded, transpose);

  /// 1. Build preliminary staggering through the x direction.

  mtk::DenseMatrix aux1(mtk::DenseMatrix::Kron(iz, iy));
  mtk::DenseMatrix dx(mtk::DenseMatrix::Kron(aux1, Dx));

  /// 2. Build preliminary staggering through the y direction.

  mtk::DenseMatrix aux2(mtk::DenseMatrix::Kron(iz, Dy));
  mtk::DenseMatrix dy(mtk::DenseMatrix::Kron(aux2, ix));

  /// 3. Build preliminary staggering through the z direction.

  mtk::DenseMatrix aux3(mtk::DenseMatrix::Kron(Dz, iy));
  mtk::DenseMatrix dz(mtk::DenseMatrix::Kron(aux3, ix));

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "Dx: " << mx << " by " << nx << std::endl;
  std::cout << "Ix: " << num_cells_x << " by " << nx  << std::endl;
  std::cout << "Dy: " << my << " by " << ny << std::endl;
  std::cout << "Iy: " << num_cells_y << " by " << ny  << std::endl;
  std::cout << "Dz: " << mz << " by " << nz << std::endl;
  std::cout << "Iz: " << num_cells_z << " by " << nz  << std::endl;
  #endif

  /// 4. Actual operator: DD_xyz = [dx dy dz].

  int total_rows{nx*ny*nz};
  int total_cols{mx*num_cells_y*num_cells_z +
                 num_cells_x*my*num_cells_z +
                 num_cells_x*num_cells_y*mz};

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "Div 3D: " << total_rows << " by " << total_cols << std::endl;
  #endif

  mtk::DenseMatrix d3d(total_rows, total_cols);

  for (auto ii = 0; ii < total_rows; ++ii) {

    for (auto jj = 0; jj < mx*num_cells_y*num_cells_z; ++jj) {
      d3d.SetValue(ii, jj, dx.GetValue(ii, jj));
    }

    int offset = mx*num_cells_y*num_cells_z;

    for(auto kk = 0; kk < num_cells_x*my*num_cells_z; ++kk) {
      d3d.SetValue(ii, kk + offset, dy.GetValue(ii, kk));
    }

    offset += num_cells_x*my*num_cells_z;

    for(auto ll = 0; ll < num_cells_x*num_cells_y*mz; ++ll) {
      d3d.SetValue(ii, ll + offset, dz.GetValue(ii, ll));
    }
  }

  divergence_ = d3d;

  return info;
}

mtk::DenseMatrix mtk::Div3D::ReturnAsDenseMatrix() const {

  return divergence_;
}
