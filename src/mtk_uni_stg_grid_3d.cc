/*!
\file mtk_uni_stg_grid_3d.cc

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

#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>

#include "mtk_tools.h"
#include "mtk_uni_stg_grid_3d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::UniStgGrid3D &in) {

  stream << '[' << in.west_bndy_ << ':' << in.num_cells_x_ << ':' <<
  in.east_bndy_ << "] x ";

  stream << '[' << in.south_bndy_ << ':' << in.num_cells_y_ << ':' <<
  in.north_bndy_ << "] x ";

  stream << '[' << in.bottom_bndy_ << ':' << in.num_cells_z_ << ':' <<
  in.top_bndy_ << "] = " << std::endl << std::endl;

  /// 1. Print spatial coordinates.

  stream << "x:";
  for (auto const &cc: in.discrete_domain_x_) {
    stream << std::setw(10) << cc;
  }
  stream << std::endl;

  stream << "y:";
  for (auto const &cc: in.discrete_domain_y_) {
    stream << std::setw(10) << cc;
  }
  stream << std::endl;

  stream << "z:";
  for (auto const &cc: in.discrete_domain_z_) {
    stream << std::setw(10) << cc;
  }
  stream << std::endl;

  /// 2. Print scalar field.

  if (in.nature_ == mtk::SCALAR) {
    stream << "u(x,y,z):" << std::endl;
    if (in.discrete_field_.size() > 0) {

    }
  } else {
    stream << "p(x,y,z):" << std::endl;
    stream << "q(x,y.z):" << std::endl;
    if (in.discrete_field_.size() > 0) {

    }
  }
  return stream;
}
}

mtk::UniStgGrid3D::UniStgGrid3D():
    discrete_domain_x_(),
    discrete_domain_y_(),
    discrete_domain_z_(),
    discrete_field_(),
    nature_(),
    west_bndy_(),
    east_bndy_(),
    num_cells_x_(),
    delta_x_(),
    south_bndy_(),
    north_bndy_(),
    num_cells_y_(),
    delta_y_(),
    bottom_bndy_(),
    top_bndy_(),
    num_cells_z_(),
    delta_z_() {}

mtk::UniStgGrid3D::UniStgGrid3D(const UniStgGrid3D &grid):
    nature_(grid.nature_),
    west_bndy_(grid.west_bndy_),
    east_bndy_(grid.east_bndy_),
    num_cells_x_(grid.num_cells_x_),
    delta_x_(grid.delta_x_),
    south_bndy_(grid.south_bndy_),
    north_bndy_(grid.north_bndy_),
    num_cells_y_(grid.num_cells_y_),
    delta_y_(grid.delta_y_),
    bottom_bndy_(grid.bottom_bndy_),
    top_bndy_(grid.top_bndy_),
    num_cells_z_(grid.num_cells_z_),
    delta_z_(grid.delta_z_) {

    std::copy(grid.discrete_domain_x_.begin(),
              grid.discrete_domain_x_.begin() + grid.discrete_domain_x_.size(),
              discrete_domain_x_.begin());

    std::copy(grid.discrete_domain_y_.begin(),
              grid.discrete_domain_y_.begin() + grid.discrete_domain_y_.size(),
              discrete_domain_y_.begin());

    std::copy(grid.discrete_domain_z_.begin(),
              grid.discrete_domain_z_.begin() + grid.discrete_domain_z_.size(),
              discrete_domain_z_.begin());

    std::copy(grid.discrete_field_.begin(),
              grid.discrete_field_.begin() + grid.discrete_field_.size(),
              discrete_field_.begin());
}

mtk::UniStgGrid3D::UniStgGrid3D(const Real &west_bndy,
                                const Real &east_bndy,
                                const int &num_cells_x,
                                const Real &south_bndy,
                                const Real &north_bndy,
                                const int &num_cells_y,
                                const Real &bottom_bndy,
                                const Real &top_bndy,
                                const int &num_cells_z,
                                const mtk::FieldNature &nature) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(west_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy <= west_bndy, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_x < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(south_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy <= south_bndy,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_y < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(bottom_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(top_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(top_bndy <= bottom_bndy,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_z < 0, __FILE__, __LINE__, __func__);
  #endif

  nature_ = nature;

  west_bndy_ = west_bndy;
  east_bndy_ = east_bndy;
  num_cells_x_ = num_cells_x;

  south_bndy_ = south_bndy;
  north_bndy_ = north_bndy;
  num_cells_y_ = num_cells_y;

  bottom_bndy_ = bottom_bndy;
  top_bndy_ = top_bndy;
  num_cells_z_ = num_cells_z;

  delta_x_ = (east_bndy_ - west_bndy_)/((mtk::Real) num_cells_x);
  delta_y_ = (north_bndy_ - south_bndy_)/((mtk::Real) num_cells_y);
  delta_z_ = (top_bndy_ - bottom_bndy_)/((mtk::Real) num_cells_z);
}

mtk::UniStgGrid3D::~UniStgGrid3D() {}

mtk::FieldNature mtk::UniStgGrid3D::nature() const {

  return nature_;
}

mtk::Real mtk::UniStgGrid3D::west_bndy() const {

  return west_bndy_;
}

mtk::Real mtk::UniStgGrid3D::east_bndy() const {

  return east_bndy_;
}

int mtk::UniStgGrid3D::num_cells_x() const {

  return num_cells_x_;
}

mtk::Real mtk::UniStgGrid3D::delta_x() const {

  return delta_x_;
}

const mtk::Real* mtk::UniStgGrid3D::discrete_domain_x() const {

  return discrete_domain_x_.data();
}

mtk::Real mtk::UniStgGrid3D::south_bndy() const {

  return south_bndy_;
}

mtk::Real mtk::UniStgGrid3D::north_bndy() const {

  return north_bndy_;
}

int mtk::UniStgGrid3D::num_cells_y() const {

  return num_cells_y_;
}

mtk::Real mtk::UniStgGrid3D::delta_y() const {

  return delta_y_;
}

const mtk::Real* mtk::UniStgGrid3D::discrete_domain_y() const {

  return discrete_domain_y_.data();
}

mtk::Real mtk::UniStgGrid3D::bottom_bndy() const {

  return bottom_bndy_;
}

mtk::Real mtk::UniStgGrid3D::top_bndy() const {

  return top_bndy_;
}

int mtk::UniStgGrid3D::num_cells_z() const {

  return num_cells_z_;
}

mtk::Real mtk::UniStgGrid3D::delta_z() const {

  return delta_z_;
}

const mtk::Real* mtk::UniStgGrid3D::discrete_domain_z() const {

  return discrete_domain_z_.data();
}

mtk::Real* mtk::UniStgGrid3D::discrete_field() {

  return discrete_field_.data();
}

bool mtk::UniStgGrid3D::Bound() const {

  return discrete_field_.size() != 0;
}

int mtk::UniStgGrid3D::Size() const {

  return discrete_field_.size();
}

void mtk::UniStgGrid3D::BindScalarField(
    mtk::Real (*ScalarField)(const mtk::Real &xx,
                             const mtk::Real &yy,
                             const mtk::Real &zz)) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nature_ != mtk::SCALAR, __FILE__, __LINE__, __func__);
  #endif

  /// 1. Create collection of spatial coordinates for \f$ x \f$.

  discrete_domain_x_.reserve(num_cells_x_ + 2);

  discrete_domain_x_.push_back(west_bndy_);
  #ifdef MTK_PRECISION_DOUBLE
  auto first_center = west_bndy_ + delta_x_/2.0;
  #else
  auto first_center = west_bndy_ + delta_x_/2.0f;
  #endif
  discrete_domain_x_.push_back(first_center);
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_domain_x_.push_back(first_center + ii*delta_x_);
  }
  discrete_domain_x_.push_back(east_bndy_);

  /// 2. Create collection of spatial coordinates for \f$ y \f$.

  discrete_domain_y_.reserve(num_cells_y_ + 2);

  discrete_domain_y_.push_back(south_bndy_);
  #ifdef MTK_PRECISION_DOUBLE
  first_center = south_bndy_ + delta_x_/2.0;
  #else
  first_center = south_bndy_ + delta_x_/2.0f;
  #endif
  discrete_domain_y_.push_back(first_center);
  for (auto ii = 1; ii < num_cells_y_; ++ii) {
    discrete_domain_y_.push_back(first_center + ii*delta_y_);
  }
  discrete_domain_y_.push_back(north_bndy_);

  /// 3. Create collection of spatial coordinates for \f$ z \f$.

  discrete_domain_z_.reserve(num_cells_z_ + 2);

  discrete_domain_z_.push_back(bottom_bndy_);
  first_center = bottom_bndy_ + delta_z_/mtk::kTwo;
  discrete_domain_z_.push_back(first_center);
  for (auto ii = 1; ii < num_cells_z_; ++ii) {
    discrete_domain_z_.push_back(first_center + ii*delta_z_);
  }
  discrete_domain_z_.push_back(top_bndy_);

  /// 4. Create collection of field samples.

  int aux{(num_cells_x_ + 2)*(num_cells_y_ + 2)*(num_cells_z_ + 2)};

  discrete_field_.reserve(aux);

  for (int kk = 0; kk < num_cells_z_ + 2; ++kk) {
    for (int ii = 0; ii < num_cells_y_ + 2; ++ii) {
      for (int jj = 0; jj < num_cells_x_ + 2; ++jj) {
        #if MTK_VERBOSE_LEVEL > 6
        std::cout << "At z = " << discrete_domain_z_[kk] << ": Pushing value"
          " for x = " << discrete_domain_x_[jj] << " y = " <<
          discrete_domain_y_[ii] << std::endl;
        #endif
        discrete_field_.push_back(ScalarField(discrete_domain_x_[jj],
                                              discrete_domain_y_[ii],
                                              discrete_domain_z_[kk]));
      }
    }
  }
}

void mtk::UniStgGrid3D::BindVectorFieldPComponent(
  mtk::Real (*VectorField)(const mtk::Real &xx,
                           const mtk::Real &yy,
                           const mtk::Real &zz)) {

}

void mtk::UniStgGrid3D::BindVectorFieldQComponent(
  mtk::Real (*VectorField)(const mtk::Real &xx,
                           const mtk::Real &yy,
                           const mtk::Real &zz)) {

}

void mtk::UniStgGrid3D::BindVectorFieldRComponent(
  mtk::Real (*VectorField)(const mtk::Real &xx,
                           const mtk::Real &yy,
                           const mtk::Real &zz)) {

}

void mtk::UniStgGrid3D::BindVectorField(
  mtk::Real (*VectorFieldPComponent)(const mtk::Real &xx,
                                     const mtk::Real &yy,
                                     const mtk::Real &zz),
  mtk::Real (*VectorFieldQComponent)(const mtk::Real &xx,
                                     const mtk::Real &yy,
                                     const mtk::Real &zz),
  mtk::Real (*VectorFieldRComponent)(const mtk::Real &xx,
                                     const mtk::Real &yy,
                                     const mtk::Real &zz)) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nature_ != mtk::VECTOR, __FILE__, __LINE__, __func__);
  #endif

  BindVectorFieldPComponent(VectorFieldPComponent);
  BindVectorFieldQComponent(VectorFieldQComponent);
}

bool mtk::UniStgGrid3D::WriteToFile(std::string filename,
                                    std::string space_name_x,
                                    std::string space_name_y,
                                    std::string space_name_z,
                                    std::string field_name) const {

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open(filename);

  if (!output_dat_file.is_open()) {
    return false;
  }

  if (nature_ == mtk::SCALAR) {
    output_dat_file << "# " << space_name_x <<  ' ' << space_name_y << ' ' <<
      space_name_z << ' ' << field_name << std::endl;

  int idx{};
  for (int kk = 0; kk < num_cells_z_ + 2; ++kk) {
    for (int ii = 0; ii < num_cells_y_ + 2; ++ii) {
      for (int jj = 0; jj < num_cells_x_ + 2; ++jj) {
        output_dat_file << discrete_domain_x_[jj] << ' ' <<
          discrete_domain_y_[ii] << ' ' << discrete_domain_z_[kk] << ' ' <<
          discrete_field_[idx] << std::endl;
        idx++;
      }
    }
  }

  } else {
    output_dat_file << "# " << space_name_x <<  ' ' << space_name_y << ' ' <<
      space_name_z << ' ' << field_name << std::endl;

  }

  output_dat_file.close();

  return true;
}
