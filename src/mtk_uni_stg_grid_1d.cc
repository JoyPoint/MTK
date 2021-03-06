/*!
\file mtk_uni_stg_grid_1d.cc

\brief Definition of an 1D uniform staggered grid.

Definition of an 1D uniform staggered grid.

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

#include <iostream>
#include <iomanip>
#include <fstream>

#include "mtk_foundations.h"
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

  /// 2. Print scalar field.

  if (in.field_nature_ == mtk::FieldNature::SCALAR) {
    stream << "u:";
  }
  else {
    stream << "v:";
  }
  for (unsigned int ii = 0; ii < in.discrete_field_.size(); ++ii) {
    stream << std::setw(10) << in.discrete_field_[ii];
  }

  stream << std::endl;

  return stream;
}
}

mtk::UniStgGrid1D& mtk::UniStgGrid1D::operator =(const mtk::UniStgGrid1D &in) {

  if(this == &in) {
    return *this;
  }

  field_nature_ = in.field_nature_;
  west_bndy_x_ = in.west_bndy_x_;
  east_bndy_x_ = in.east_bndy_x_;
  num_cells_x_ = in.num_cells_x_;
  delta_x_ = in.delta_x_;

  discrete_domain_x_.clear();

  std::copy(in.discrete_domain_x_.begin(),
            in.discrete_domain_x_.end(),
            discrete_domain_x_.begin());

  discrete_field_.clear();

  std::copy(in.discrete_field_.begin(),
            in.discrete_field_.end(),
            discrete_field_.begin());

  return *this;
}

mtk::UniStgGrid1D::UniStgGrid1D():
    field_nature_(),
    discrete_domain_x_(),
    discrete_field_(),
    west_bndy_x_(),
    east_bndy_x_(),
    num_cells_x_(),
    delta_x_() {}

mtk::UniStgGrid1D::UniStgGrid1D(const UniStgGrid1D &in):
    field_nature_(in.field_nature_),
    west_bndy_x_(in.west_bndy_x_),
    east_bndy_x_(in.east_bndy_x_),
    num_cells_x_(in.num_cells_x_),
    delta_x_(in.delta_x_) {

    std::copy(in.discrete_domain_x_.begin(),
              in.discrete_domain_x_.begin() + in.discrete_domain_x_.size(),
              discrete_domain_x_.begin());

    std::copy(in.discrete_field_.begin(),
              in.discrete_field_.begin() + in.discrete_field_.size(),
              discrete_field_.begin());
}

mtk::UniStgGrid1D::UniStgGrid1D(const Real &west_bndy_x,
                                const Real &east_bndy_x,
                                const int &num_cells_x,
                                const mtk::FieldNature &field_nature) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(west_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy_x <= west_bndy_x, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_x < 0, __FILE__, __LINE__, __func__);
  #endif

  field_nature_ = field_nature;
  west_bndy_x_ = west_bndy_x;
  east_bndy_x_ = east_bndy_x;
  num_cells_x_ = num_cells_x;

  delta_x_ = (east_bndy_x - west_bndy_x)/((mtk::Real) num_cells_x);
}

mtk::UniStgGrid1D::~UniStgGrid1D() {}

mtk::Real mtk::UniStgGrid1D::west_bndy_x() const {

  return west_bndy_x_;
}

mtk::Real mtk::UniStgGrid1D::east_bndy_x() const {

  return east_bndy_x_;
}

mtk::Real mtk::UniStgGrid1D::delta_x() const {

  return delta_x_;
}

const mtk::Real *mtk::UniStgGrid1D::discrete_domain_x() const {

  return discrete_domain_x_.data();
}

mtk::Real *mtk::UniStgGrid1D::discrete_field() {

  return discrete_field_.data();
}

int mtk::UniStgGrid1D::num_cells_x() const {

  return num_cells_x_;
}

mtk::FieldNature mtk::UniStgGrid1D::field_nature() const {

  return field_nature_;
}

void mtk::UniStgGrid1D::GenerateDiscreteDomainX() {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(discrete_domain_x_.size() != 0,
                      __FILE__, __LINE__, __func__);
  #endif

  if (field_nature_ == mtk::FieldNature::SCALAR) {

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

  } else {

    discrete_domain_x_.reserve(num_cells_x_ + 1);

    discrete_domain_x_.push_back(west_bndy_x_);
    for (auto ii = 1; ii < num_cells_x_; ++ii) {
      discrete_domain_x_.push_back(west_bndy_x_ + ii*delta_x_);
    }
    discrete_domain_x_.push_back(east_bndy_x_);
  }
}

void mtk::UniStgGrid1D::ReserveDiscreteField() {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(discrete_field_.size() != 0,
                      __FILE__, __LINE__, __func__);
  #endif

  if (field_nature_ == mtk::FieldNature::SCALAR) {
    discrete_field_.reserve(num_cells_x_ + 2);
  } else {
    discrete_field_.reserve(num_cells_x_ + 1);
  }
}

void mtk::UniStgGrid1D::BindScalarField(
    mtk::Real (*ScalarField)(const mtk::Real &xx,
                             const std::vector<mtk::Real> &pp),
    const std::vector<Real> &parameters) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(field_nature_ == mtk::FieldNature::VECTOR,
                      __FILE__, __LINE__, __func__);
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

  std::vector<mtk::Real> aux(parameters);

  discrete_field_.reserve(num_cells_x_ + 2);

  discrete_field_.push_back(ScalarField(west_bndy_x_, aux));

  discrete_field_.push_back(ScalarField(first_center, aux));
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_field_.push_back(ScalarField(first_center + ii*delta_x_,
                                          aux));
  }
  discrete_field_.push_back(ScalarField(east_bndy_x_, aux));
}

void mtk::UniStgGrid1D::BindScalarField(const std::vector<Real> &samples) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(field_nature_ == mtk::FieldNature::VECTOR,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(samples.size() != num_cells_x_,
                      __FILE__, __LINE__, __func__);
  #endif

  discrete_field_ = samples;
}

void mtk::UniStgGrid1D::BindVectorField(
    mtk::Real (*VectorField)(const mtk::Real &xx,
                             const std::vector<mtk::Real> &pp),
    const std::vector<Real> &parameters) {

  #ifdef MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(field_nature_ == mtk::FieldNature::SCALAR,
                      __FILE__, __LINE__, __func__);
  #endif

  /// 1. Create collection of spatial coordinates.

  discrete_domain_x_.reserve(num_cells_x_ + 1);

  discrete_domain_x_.push_back(west_bndy_x_);
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_domain_x_.push_back(west_bndy_x_ + ii*delta_x_);
  }
  discrete_domain_x_.push_back(east_bndy_x_);

  /// 2. Create collection of field samples.

  std::vector<mtk::Real> aux(parameters);

  discrete_field_.reserve(num_cells_x_ + 1);

  discrete_field_.push_back(VectorField(west_bndy_x_, aux));
  for (auto ii = 1; ii < num_cells_x_; ++ii) {
    discrete_field_.push_back(VectorField(
      west_bndy_x_ + ii*delta_x_, aux));
  }
  discrete_field_.push_back(VectorField(east_bndy_x_, aux));
}

bool mtk::UniStgGrid1D::WriteToFile(std::string filename,
                                    std::string space_name,
                                    std::string field_name) const {

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open(filename);

  if (!output_dat_file.is_open()) {
    return false;
  }

  output_dat_file << "# " << space_name <<  ' ' << field_name << std::endl;
  for (unsigned int ii = 0; ii < discrete_domain_x_.size(); ++ii) {
    output_dat_file << discrete_domain_x_[ii] << ' ' << discrete_field_[ii] <<
      std::endl;
  }

  output_dat_file.close();

  return true;
}
