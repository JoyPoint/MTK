/*!
\file mtk_div_2d.cc

\brief Implements the class Div2D.

This class implements a 2D divergence matrix operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm.

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
#include "mtk_enums.h"
#include "mtk_uni_stg_grid_1d.h"
#include "mtk_div_1d.h"
#include "mtk_div_2d.h"

mtk::Div2D::Div2D():
  order_accuracy_(),
  mimetic_threshold_() {}

mtk::Div2D::Div2D(const Div2D &div):
  order_accuracy_(div.order_accuracy_),
  mimetic_threshold_(div.mimetic_threshold_) {}

mtk::Div2D::~Div2D() {}

bool mtk::Div2D::ConstructDiv2D(const mtk::UniStgGrid2D &grid,
                                int order_accuracy,
                                mtk::Real mimetic_threshold) {

  int num_cells_x = grid.num_cells_x();
  int num_cells_y = grid.num_cells_y();

  int mx = num_cells_x + 2;  // Dx vertical dimension.
  int nx = num_cells_x + 1;  // Dx horizontal dimension.
  int my = num_cells_y + 2;  // Dy vertical dimension.
  int ny = num_cells_y + 1;  // Dy horizontal dimension.

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

  mtk::UniStgGrid1D grid_x(west, east, num_cells_x, mtk::FieldNature::VECTOR);
  mtk::UniStgGrid1D grid_y(south, north, num_cells_y, mtk::FieldNature::VECTOR);

  mtk::DenseMatrix dx(div.ReturnAsDenseMatrix(grid_x));
  mtk::DenseMatrix dy(div.ReturnAsDenseMatrix(grid_y));

  bool padded{true};
  bool transpose{false};

  mtk::DenseMatrix ix(num_cells_x, padded, transpose);
  mtk::DenseMatrix iy(num_cells_y, padded, transpose);

  mtk::DenseMatrix dxy(mtk::DenseMatrix::Kron(iy, dx));
  mtk::DenseMatrix dyx(mtk::DenseMatrix::Kron(dy, ix));

  #if MTK_VERBOSE_LEVEL > 2
  std::cout << "Dx: " << mx << " by " << nx << std::endl;
  std::cout << "Iy : " << num_cells_y<< " by " << ny  << std::endl;
  std::cout << "Dy: " << my << " by " << ny << std::endl;
  std::cout << "Ix : " << num_cells_x<< " by " << nx  << std::endl;
  std::cout << "Div 2D: " << mx*num_cells_y + my*num_cells_x << " by " <<
    nx*ny <<std::endl;
  #endif

  mtk::DenseMatrix d2d(mx*my, nx*num_cells_y + ny*num_cells_x);

  for (auto ii = 0; ii < mx*my; ii++) {
    for (auto jj = 0; jj < nx*num_cells_y; jj++) {
      d2d.SetValue(ii, jj, dxy.GetValue(ii,jj));
    }
    for(auto kk = 0; kk<ny*num_cells_x; kk++) {
      d2d.SetValue(ii, kk + nx*num_cells_y, dyx.GetValue(ii, kk));
    }
  }

  divergence_ = d2d;

  return info;
}

mtk::DenseMatrix mtk::Div2D::ReturnAsDenseMatrix() const {

  return divergence_;
}
