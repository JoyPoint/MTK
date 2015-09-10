/*!
\file mtk_uni_stg_grid_1d.cc

\brief Implementation of an 1D uniform staggered grid.

Implementation of an 1D uniform staggered grid.

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
#include <fstream>

#include "mtk_roots.h"
#include "mtk_enums.h"
#include "mtk_tools.h"

#include "mtk_uni_stg_grid_1d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::UniStgGrid1D &in) {

  stream << '[' << in.west_bndy_x_ << ':' << in.num_cells_x_ << ':' <<
  in.east_bndy_x_ << "] = " << std::endl << std::endl;

  /// 1. Print spatial coordinates.

  stream << "x:";
  for (unsigned int ii = 0; ii < in.discrete_domain_x_.size(); ++ii) {
    stream << std::setw(10) << in.discrete_domain_x_[ii];
  }
  stream << std::endl;

  if (in.nature_ == mtk::SCALAR) {
    stream << "u:";
  }
  else {
    stream << "v:";
  }
  for (unsigned int ii = 0; ii < in.discrete_field_u_.size(); ++ii) {
    stream << std::setw(10) << in.discrete_field_u_[ii];
  }

  stream << std::endl;

  return stream;
}
}

mtk::UniStgGrid1D::UniStgGrid1D():
    nature_(),
    discrete_domain_x_(),
    discrete_field_u_(),
    west_bndy_x_(),
    east_bndy_x_(),
    num_cells_x_(),
    delta_x_() {}

mtk::UniStgGrid1D::UniStgGrid1D(const UniStgGrid1D &grid):
    nature_(grid.nature_),
    west_bndy_x_(grid.west_bndy_x_),
    east_bndy_x_(grid.east_bndy_x_),
    num_cells_x_(grid.num_cells_x_),
    delta_x_(grid.delta_x_) {

    std::copy(grid.discrete_domain_x_.begin(),
              grid.discrete_domain_x_.begin() + grid.discrete_domain_x_.size(),
              discrete_domain_x_.begin());

    std::copy(grid.discrete_field_u_.begin(),
              grid.discrete_field_u_.begin() + grid.discrete_field_u_.size(),
              discrete_field_u_.begin());
}

mtk::UniStgGrid1D::UniStgGrid1D(const Real &west_bndy_x,
                                const Real &east_bndy_x,
                                const int &num_cells_x,
                                const mtk::FieldNature &nature) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x <= west_bndy_x, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_x < 0, __FILE__, __LINE__, __func__);
  #endif

  nature_ = nature;
  west_bndy_x_ = west_bndy_x;
  east_bndy_x_ = east_bndy_x;
  num_cells_x_ = num_cells_x;

  delta_x_ = (east_bndy_x - west_bndy_x)/((mtk::Real) num_cells_x);
}

mtk::UniStgGrid1D::~UniStgGrid1D() {}

mtk::Real mtk::UniStgGrid1D::delta_x() const {

  return delta_x_;
}

mtk::Real *mtk::UniStgGrid1D::discrete_domain_x() {

  return discrete_domain_x_.data();
}

mtk::Real *mtk::UniStgGrid1D::discrete_field_u() {

  return discrete_field_u_.data();
}

int mtk::UniStgGrid1D::num_cells_x() const {

  return num_cells_x_;
}

void mtk::UniStgGrid1D::BindScalarField(
    mtk::Real (*ScalarField)(mtk::Real xx)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(nature_ == mtk::VECTOR, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Create collection of spatial coordinates.

  discrete_domain_x_.reserve(num_cells_x_ + 2);

  discrete_domain_x_.push_back(west_bndy_x_);
  #ifdef MTK_PRECISION_DOUBLE
  auto first_center = west_bndy_x_ + delta_x_/2.0;
  #else
  auto first_center = west_bndy_x_ + delta_x_/2.0f;
  #endif
  discrete_domain_x_.push_back(first_center);
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_domain_x_.push_back(first_center + ii*delta_x_);
  }
  discrete_domain_x_.push_back(east_bndy_x_);

  /// 2. Create collection of field samples.

  discrete_field_u_.reserve(num_cells_x_ + 2);

  discrete_field_u_.push_back(ScalarField(west_bndy_x_));

  discrete_field_u_.push_back(ScalarField(first_center));
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_field_u_.push_back(ScalarField(first_center + ii*delta_x_));
  }
  discrete_field_u_.push_back(ScalarField(east_bndy_x_));
}

void mtk::UniStgGrid1D::BindVectorField(
    mtk::Real (*VectorField)(mtk::Real xx)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(nature_ == mtk::SCALAR, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Create collection of spatial coordinates.

  discrete_domain_x_.reserve(num_cells_x_ + 1);

  discrete_domain_x_.push_back(west_bndy_x_);
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_domain_x_.push_back(west_bndy_x_ + ii*delta_x_);
  }
  discrete_domain_x_.push_back(east_bndy_x_);

  /// 2. Create collection of field samples.

  discrete_field_u_.reserve(num_cells_x_ + 1);

  discrete_field_u_.push_back(VectorField(west_bndy_x_));
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_field_u_.push_back(VectorField(west_bndy_x_ + ii*delta_x_));
  }
  discrete_field_u_.push_back(VectorField(east_bndy_x_));
}

bool mtk::UniStgGrid1D::WriteToFile(std::string filename,
                                    std::string space_name,
                                    std::string field_name) {

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open(filename);

  if (!output_dat_file.is_open()) {
    return false;
  }

  output_dat_file << "# " << space_name <<  ' ' << field_name << std::endl;
  for (unsigned int ii = 0; ii < discrete_domain_x_.size(); ++ii) {
    output_dat_file << discrete_domain_x_[ii] << ' ' << discrete_field_u_[ii] <<
      std::endl;
  }

  output_dat_file.close();

  return true;
}
