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
#include "mtk_uni_stg_grid_2d.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::UniStgGrid2D &in) {

  stream << '[' << in.west_bndy_ << ':' << in.num_cells_x_ << ':' <<
  in.east_bndy_ << "] x ";

  stream << '[' << in.south_bndy_ << ':' << in.num_cells_y_ << ':' <<
  in.north_bndy_ << "] = " << std::endl << std::endl;

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
    stream << "u:" << std::endl;
    if (in.discrete_field_.size() > 0) {
      for (int ii = 0; ii < in.num_cells_x_ + 2; ++ii) {
        for (int jj = 0; jj < in.num_cells_y_ + 2; ++jj) {
          stream << std::setw(10) << in.discrete_field_[ii*in.num_cells_y_ + jj];
        }
        stream << std::endl;
      }
    }
  } else {

    int mm{in.num_cells_x_};
    int nn{in.num_cells_y_};
    int p_offset{nn*(mm + 1) - 1};

    stream << "p(x,y):" << std::endl;
    for (int ii = 0; ii < nn; ++ii) {
      for (int jj = 0; jj < mm + 1; ++jj) {
        stream << std::setw(10) << in.discrete_field_[ii*(mm + 1) + jj];
      }
      stream << std::endl;
    }
    stream << std::endl;

    stream << "q(x,y):" << std::endl;
    for (int ii = 0; ii < nn + 1; ++ii) {
      for (int jj = 0; jj < mm; ++jj) {
        stream << std::setw(10) <<
          in.discrete_field_[p_offset + ii*mm + jj];
      }
      stream << std::endl;
    }
    stream << std::endl;
  }

  return stream;
}
}

mtk::UniStgGrid2D::UniStgGrid2D():
    discrete_domain_x_(),
    discrete_domain_y_(),
    discrete_field_(),
    nature_(),
    west_bndy_(),
    east_bndy_(),
    num_cells_x_(),
    delta_x_(),
    south_bndy_(),
    north_bndy_(),
    num_cells_y_(),
    delta_y_() {}

mtk::UniStgGrid2D::UniStgGrid2D(const UniStgGrid2D &grid):
    nature_(grid.nature_),
    west_bndy_(grid.west_bndy_),
    east_bndy_(grid.east_bndy_),
    num_cells_x_(grid.num_cells_x_),
    delta_x_(grid.delta_x_),
    south_bndy_(grid.south_bndy_),
    north_bndy_(grid.north_bndy_),
    num_cells_y_(grid.num_cells_y_),
    delta_y_(grid.delta_y_) {

    std::copy(grid.discrete_domain_x_.begin(),
              grid.discrete_domain_x_.begin() + grid.discrete_domain_x_.size(),
              discrete_domain_x_.begin());

    std::copy(grid.discrete_domain_y_.begin(),
              grid.discrete_domain_y_.begin() + grid.discrete_domain_y_.size(),
              discrete_domain_y_.begin());

    std::copy(grid.discrete_field_.begin(),
              grid.discrete_field_.begin() + grid.discrete_field_.size(),
              discrete_field_.begin());
}

mtk::UniStgGrid2D::UniStgGrid2D(const Real &west_bndy,
                                const Real &east_bndy,
                                const int &num_cells_x,
                                const Real &south_bndy,
                                const Real &north_bndy,
                                const int &num_cells_y,
                                const mtk::FieldNature &nature) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(west_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east_bndy <= west_bndy, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_x < 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(south_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy < mtk::kZero, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(north_bndy <= south_bndy,
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(num_cells_y < 0, __FILE__, __LINE__, __func__);
  #endif

  nature_ = nature;

  west_bndy_ = west_bndy;
  east_bndy_ = east_bndy;
  num_cells_x_ = num_cells_x;

  south_bndy_ = south_bndy;
  north_bndy_ = north_bndy;
  num_cells_y_ = num_cells_y;

  delta_x_ = (east_bndy_ - west_bndy_)/((mtk::Real) num_cells_x);
  delta_y_ = (north_bndy_ - south_bndy_)/((mtk::Real) num_cells_y);
}

mtk::UniStgGrid2D::~UniStgGrid2D() {}

mtk::FieldNature mtk::UniStgGrid2D::nature() const {

  return nature_;
}

mtk::Real mtk::UniStgGrid2D::west_bndy() const {

  return west_bndy_;
}

mtk::Real mtk::UniStgGrid2D::east_bndy() const {

  return east_bndy_;
}

int mtk::UniStgGrid2D::num_cells_x() const {

  return num_cells_x_;
}

mtk::Real mtk::UniStgGrid2D::delta_x() const {

  return delta_x_;
}

mtk::Real mtk::UniStgGrid2D::south_bndy() const {

  return south_bndy_;
}

mtk::Real mtk::UniStgGrid2D::north_bndy() const {

  return north_bndy_;
}

int mtk::UniStgGrid2D::num_cells_y() const {

  return num_cells_y_;
}

mtk::Real mtk::UniStgGrid2D::delta_y() const {

  return delta_y_;
}

void mtk::UniStgGrid2D::BindScalarField(Real (*ScalarField)(Real xx, Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
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

  /// 3. Create collection of field samples.

  discrete_field_.reserve((num_cells_x_ + 2)*(num_cells_y_ + 2));

  for (int ii = 0; ii < num_cells_x_ + 2; ++ii) {
    for (int jj = 0; jj < num_cells_y_ + 2; ++jj) {
      discrete_field_.push_back(ScalarField(discrete_domain_x_[ii],
                                            discrete_domain_y_[jj]));
    }
  }
}

void mtk::UniStgGrid2D::BindVectorFieldPComponent(
  mtk::Real (*VectorField)(mtk::Real xx, mtk::Real yy)) {

  int mm{num_cells_x_};
  int nn{num_cells_y_};

  int total{nn*(mm + 1) + mm*(nn + 1)};

  #ifdef MTK_PRECISION_DOUBLE
  double half_delta_x{delta_x_/2.0};
  double half_delta_y{delta_y_/2.0};
  #else
  float half_delta_x{delta_x_/2.0f};
  float half_delta_y{delta_y_/2.0f};
  #endif

  /// 1. Create collection of spatial coordinates for \f$ x \f$.

  // We need every data point of the discrete domain; i.e. we need all the
  // nodes and all the centers. There are mm centers for the x direction, and
  // nn centers for the y direction. Since there is one node per center, that
  // amounts to 2*mm. If we finally consider the final boundary node, it
  // amounts to a total of 2*mm + 1 for the x direction. Analogously, for the
  // y direction, this amounts to 2*nn + 1.

  discrete_domain_x_.reserve(2*mm + 1);

  discrete_domain_x_.push_back(west_bndy_);
  for (int ii = 1; ii < (2*mm + 1); ++ii) {
    discrete_domain_x_.push_back(west_bndy_ + ii*half_delta_x);
  }

  /// 2. Create collection of spatial coordinates for \f$ y \f$.

  discrete_domain_y_.reserve(2*nn + 1);

  discrete_domain_y_.push_back(south_bndy_);
  for (int ii = 1; ii < (2*nn + 1); ++ii) {
    discrete_domain_y_.push_back(south_bndy_ + ii*half_delta_y);
  }

  /// 3. Allocate space for discrete vector field and bind \$ p \$ component.

  discrete_field_.reserve(total);

  // For each y-center.
  for (int ii = 1; ii < 2*nn + 1; ii += 2) {

    // Bind all of the x-nodes for this y-center.
    for (int jj = 0; jj < 2*mm + 1; jj += 2) {
      discrete_field_.push_back(VectorField(discrete_domain_x_[jj],
                                            discrete_domain_y_[ii]));

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Binding v at x = " << discrete_domain_x_[jj] << " y = " <<
        discrete_domain_y_[ii] << " = " <<
        VectorField(discrete_domain_x_[jj],discrete_domain_y_[ii]) << std::endl;
      #endif
    }
  }
  #if MTK_DEBUG_LEVEL > 0
  std::cout << std::endl;
  #endif
}

void mtk::UniStgGrid2D::BindVectorFieldQComponent(
  mtk::Real (*VectorField)(mtk::Real xx, mtk::Real yy)) {

  int mm{num_cells_x_};
  int nn{num_cells_y_};

  /// 3. Bind \$ q \$ component, since \$ p \$ component has already been bound.

  // For each y-node.
  for (int ii = 0; ii < 2*nn + 1; ii += 2) {

    // Bind all of the x-center for this y-node.
    for (int jj = 1; jj < 2*mm + 1; jj += 2) {
      discrete_field_.push_back(VectorField(discrete_domain_x_[jj],
                                            discrete_domain_y_[ii]));

      #if MTK_DEBUG_LEVEL > 0
      std::cout << "Binding v at x = " << discrete_domain_x_[jj] << " y = " <<
        discrete_domain_y_[ii] << " = " <<
        VectorField(discrete_domain_x_[jj],discrete_domain_y_[ii]) << std::endl;
      #endif
    }
  }
  #if MTK_DEBUG_LEVEL > 0
  std::cout << std::endl;
  #endif
}

void mtk::UniStgGrid2D::BindVectorField(
  Real (*VectorFieldPComponent)(Real xx,Real yy),
  Real (*VectorFieldQComponent)(Real xx,Real yy)) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(nature_ != mtk::VECTOR, __FILE__, __LINE__, __func__);
  #endif

  BindVectorFieldPComponent(VectorFieldPComponent);
  BindVectorFieldQComponent(VectorFieldQComponent);
}

bool mtk::UniStgGrid2D::WriteToFile(std::string filename,
                                    std::string space_name_x,
                                    std::string space_name_y,
                                    std::string field_name) {

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open(filename);

  if (!output_dat_file.is_open()) {
    return false;
  }

  if (nature_ == mtk::SCALAR) {
    output_dat_file << "# " << space_name_x <<  ' ' << space_name_y << ' ' <<
      field_name << std::endl;

    for (unsigned int ii = 0; ii < discrete_domain_x_.size(); ++ii) {
      for (unsigned int jj = 0; jj < discrete_domain_y_.size(); ++jj) {
        output_dat_file << discrete_domain_x_[ii] << ' ' <<
                           discrete_domain_y_[jj] << ' ' <<
                           discrete_field_[ii*discrete_domain_y_.size() + jj] <<
                          std::endl;
      }
      output_dat_file << std::endl;
    }
  } else {
    output_dat_file << "# " << space_name_x <<  ' ' << space_name_y << ' ' <<
      field_name << std::endl;

    output_dat_file << "# Horizontal component:" << std::endl;

    int mm{num_cells_x_};
    int nn{num_cells_y_};

    /// Write the values of the p component, with a null q component.

    // For each y-center.
    int idx{};
    for (int ii = 1; ii < 2*nn + 1; ii += 2) {
      // Bind all of the x-nodes for this y-center.
      for (int jj = 0; jj < 2*mm + 1; jj += 2) {

        output_dat_file << discrete_domain_x_[jj] << ' ' <<
          discrete_domain_y_[ii] << ' ' << discrete_field_[idx] << ' ' <<
          mtk::kZero << std::endl;

        ++idx;
      }
    }

    /// Write the values of the q component, with a null p component.
    int p_offset{nn*(mm + 1) - 1};
    idx = 0;
    output_dat_file << "# Vertical component:" << std::endl;
    // For each y-node.
    for (int ii = 0; ii < 2*nn + 1; ii += 2) {
      // Bind all of the x-center for this y-node.
      for (int jj = 1; jj < 2*mm + 1; jj += 2) {

        output_dat_file << discrete_domain_x_[jj] << ' ' <<
          discrete_domain_y_[ii] << ' ' << mtk::kZero << ' ' <<
          discrete_field_[p_offset + idx] << std::endl;

        ++idx;
      }
    }
  }

  output_dat_file.close();

  return true;
}
