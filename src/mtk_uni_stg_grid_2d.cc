/*!
\file mtk_uni_stg_grid_2d.cc

\brief Implementation of a 2D uniform staggered grid.

Implementation of a 2D uniform staggered grid.

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

#include <iostream>
#include <iomanip>

#include <algorithm>

#include "mtk_tools.h"
#include "mtk_uni_stg_grid_2d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::UniStgGrid2D &in) {

  stream << '[' << in.west_bndy_x_ << ':' << in.num_cells_x_ << ':' <<
  in.east_bndy_x_ << "] x ";

  stream << '[' << in.south_bndy_y_ << ':' << in.num_cells_y_ << ':' <<
  in.north_bndy_y_ << "] = " << std::endl << std::endl;

  /// 1. Print spatial coordinates.

  stream << "x:";
  for (unsigned int ii = 0; ii < in.discrete_domain_x_.size(); ++ii) {
    stream << std::setw(10) << in.discrete_domain_x_[ii];
  }
  stream << std::endl;

  stream << "y:";
  for (unsigned int ii = 0; ii < in.discrete_domain_y_.size(); ++ii) {
    stream << std::setw(10) << in.discrete_domain_y_[ii];
  }
  stream << std::endl;

  /// 2. Print scalar field.

  if (in.nature_ == mtk::SCALAR) {
    stream << "u:";
  }
  else {
    stream << "v:";
  }
  if (in.discrete_field_u_.size() > 0) {
    for (unsigned int ii = 0; ii < in.num_cells_x_ + 2; ++ii) {
      for (unsigned int jj = 0; jj < in.num_cells_y_ + 2; ++jj) {
        stream << std::setw(10) << in.discrete_field_u_[ii*in.num_cells_y_ + jj];
      }
      stream << std::endl;
    }
  }

  stream << std::endl;

  return stream;
}
}

mtk::UniStgGrid2D::UniStgGrid2D():
    nature_(),
    discrete_domain_x_(),
    discrete_domain_y_(),
    discrete_field_u_(),
    west_bndy_x_(),
    east_bndy_x_(),
    num_cells_x_(),
    delta_x_(),
    south_bndy_y_(),
    north_bndy_y_(),
    num_cells_y_(),
    delta_y_() {}

mtk::UniStgGrid2D::UniStgGrid2D(const UniStgGrid2D &grid):
    nature_(grid.nature_),
    west_bndy_x_(grid.west_bndy_x_),
    east_bndy_x_(grid.east_bndy_x_),
    num_cells_x_(grid.num_cells_x_),
    delta_x_(grid.delta_x_),
    south_bndy_y_(grid.south_bndy_y_),
    north_bndy_y_(grid.north_bndy_y_),
    num_cells_y_(grid.num_cells_y_),
    delta_y_(grid.delta_y_) {

    std::copy(grid.discrete_domain_x_.begin(),
              grid.discrete_domain_x_.begin() + grid.discrete_domain_x_.size(),
              discrete_domain_x_.begin());

    std::copy(grid.discrete_domain_y_.begin(),
              grid.discrete_domain_y_.begin() + grid.discrete_domain_y_.size(),
              discrete_domain_y_.begin());

    std::copy(grid.discrete_field_u_.begin(),
              grid.discrete_field_u_.begin() + grid.discrete_field_u_.size(),
              discrete_field_u_.begin());
}

mtk::UniStgGrid2D::UniStgGrid2D(const Real &west_bndy_x,
                                const Real &east_bndy_x,
                                const int &num_cells_x,
                                const Real &south_bndy_y,
                                const Real &north_bndy_y,
                                const int &num_cells_y,
                                const mtk::FieldNature &nature) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x <= west_bndy_x, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_x < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(south_bndy_y < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy_y < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy_y <= south_bndy_y,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_y < 0, __FILE__, __LINE__, __func__);
  #endif

  nature_ = nature;

  west_bndy_x_ = west_bndy_x;
  east_bndy_x_ = east_bndy_x;
  num_cells_x_ = num_cells_x;

  south_bndy_y_ = south_bndy_y;
  north_bndy_y_ = north_bndy_y;
  num_cells_y_ = num_cells_y;

  delta_x_ = (east_bndy_x - west_bndy_x)/((mtk::Real) num_cells_x);
  delta_y_ = (north_bndy_y - south_bndy_y)/((mtk::Real) num_cells_y);
}

mtk::UniStgGrid2D::~UniStgGrid2D() {}

mtk::Real mtk::UniStgGrid2D::west_bndy_x() const {

  return west_bndy_x_;
}

mtk::Real mtk::UniStgGrid2D::east_bndy_x() const {

  return east_bndy_x_;
}

mtk::Real mtk::UniStgGrid2D::south_bndy_y() const {

  return south_bndy_y_;
}

mtk::Real mtk::UniStgGrid2D::north_bndy_y() const {

  return north_bndy_y_;
}
