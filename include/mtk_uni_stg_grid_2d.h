/*!
\file mtk_uni_stg_grid_2d.h

\brief Definition of an 2D uniform staggered grid.

Definition of an 1D uniform staggered grid.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\todo Create overloaded binding routines that read data from files.
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

#ifndef MTK_INCLUDE_UNI_STG_GRID_2D_H_
#define MTK_INCLUDE_UNI_STG_GRID_2D_H_

#include <vector>
#include <string>

#include "mtk_roots.h"
#include "mtk_enums.h"
#include "mtk_dense_matrix.h"

namespace mtk {

/*!
\class UniStgGrid2D

\ingroup c06-grids

\brief Uniform 2D Staggered Grid.

Uniform 2D Staggered Grid.
*/
class UniStgGrid2D {
 public:
  /// \brief Prints the grid as a tuple of arrays.
  friend std::ostream& operator <<(std::ostream& stream, UniStgGrid2D &in);

  /// \brief Default constructor.
  UniStgGrid2D();

  /*!
  \brief Copy constructor.

  \param [in] grid Given grid.
  */
  UniStgGrid2D(const UniStgGrid2D &grid);

  /*!
  \brief Construct a grid based on spatial discretization parameters.

  \param[in] west_bndy_x Coordinate for the west boundary.
  \param[in] east_bndy_x Coordinate for the east boundary.
  \param[in] num_cells_x Number of cells of the required grid.
  \param[in] south_bndy_y Coordinate for the west boundary.
  \param[in] north_bndy_y Coordinate for the east boundary.
  \param[in] num_cells_y Number of cells of the required grid.
  \param[in] nature Nature of the discrete field to hold.

  \sa mtk::FieldNature
  */
  UniStgGrid2D(const Real &west_bndy_x,
               const Real &east_bndy_x,
               const int &num_cells_x,
               const Real &south_bndy_y,
               const Real &north_bndy_y,
               const int &num_cells_y,
               const mtk::FieldNature &nature = mtk::SCALAR);

  /// \brief Destructor.
  ~UniStgGrid2D();

  /*!
  \brief Provides access to the grid spatial data.

  \return Pointer to the spatial data.
  */
  Real *discrete_domain_x();

  /*!
  \brief Provides access to the grid spatial data.

  \return Pointer to the spatial data.
  */
  Real *discrete_domain_y();

  /*!
  \brief Provides access to the grid field data.

  \return Pointer to the field data.
  */
  Real *discrete_field();

  /*!
  \brief Physical nature of the data bound to the grid.

  \return Value of an enumeration.

  \sa mtk::FieldNature
  */
  FieldNature nature() const;

  /*!
  \brief Provides access to west boundary spatial coordinate.

  \return West boundary spatial coordinate.
  */
  Real west_bndy() const;

  /*!
  \brief Provides access to east boundary spatial coordinate.

  \return East boundary spatial coordinate.
  */
  Real east_bndy() const;

  /*!
  \brief Provides access to the number of cells of the grid.

  \return Number of cells of the grid.
  */
  int num_cells_x() const;

  /*!
  \brief Provides access to the computed \$ \Delta x \$.

  \return Computed \$ \Delta x \$.
  */
  Real delta_x() const;

  /*!
  \brief Provides access to south boundary spatial coordinate.

  \return South boundary spatial coordinate.
  */
  Real south_bndy() const;

  /*!
  \brief Provides access to north boundary spatial coordinate.

  \return North boundary spatial coordinate.
  */
  Real north_bndy() const;

  /*!
  \brief Provides access to the number of cells of the grid.

  \return Number of cells of the grid.
  */
  int num_cells_y() const;

  /*!
  \brief Provides access to the computed \$ \Delta y \$.

  \return Computed \$ \Delta y \$.
  */
  Real delta_y() const;

  /*!
  \brief Binds a given scalar field to the grid.

  \param[in] ScalarField Pointer to the function implementing the scalar field.
  */
  void BindScalarField(Real (*ScalarField)(Real xx, Real yy));

  /*!
  \brief Binds a given vector field to the grid.

  We assume the field to be of the form:

  \f[
    \mathbf{v}(x) = p(x, y)\hat{\mathbf{i}} + q(x, y)\hat{\mathbf{j}}
  \f]

  \param[in] VectorFieldPComponent Pointer to the function implementing the \$
  p \$ component of the vector field.
  \param[in] VectorFieldPComponent Pointer to the function implementing the \$
  q \$ component of the vector field.
  */
  void BindVectorField(Real (*VectorFieldPComponent)(Real xx,Real yy),
                       Real (*VectorFieldQComponent)(Real xx,Real yy));

  /*!
  \brief Writes grid to a file compatible with Gnuplot 4.6.

  \param[in] filename Name of the output file.
  \param[in] space_name_x Name for the first column of the (spatial) data.
  \param[in] space_name_y Name for the second column of the (spatial) data.
  \param[in] field_name Name for the second column of the (physical field) data.

  \return Success of the file writing process.

  \sa http://www.gnuplot.info/
  */
  bool WriteToFile(std::string filename,
                   std::string space_name_x,
                   std::string space_name_y,
                   std::string field_name);

 private:
  /*!
  \brief Binds a given component of a vector field to the grid.

  We assume the field to be of the form:

  \f[
    \mathbf{v}(x) = p(x, y)\hat{\mathbf{i}} + q(x, y)\hat{\mathbf{j}}
  \f]

  \param[in] BindVectorFieldPComponent Pointer to the function implementing the
  \$ p \$ component of the vector field.
  */
  void BindVectorFieldPComponent(
    Real (*VectorFieldPComponent)(Real xx, Real yy));

  /*!
  \brief Binds a given component of a vector field to the grid.

  We assume the field to be of the form:

  \f[
    \mathbf{v}(x) = p(x, y)\hat{\mathbf{i}} + q(x, y)\hat{\mathbf{j}}
  \f]

  \param[in] BindVectorFieldQComponent Pointer to the function implementing the
  \$ q \$ component of the vector field.
  */
  void BindVectorFieldQComponent(
    Real (*VectorFieldQComponent)(Real xx, Real yy));


  std::vector<Real> discrete_domain_x_; ///< Array of spatial data.
  std::vector<Real> discrete_domain_y_; ///< Array of spatial data.
  std::vector<Real> discrete_field_;    ///< Array of field's data.

  FieldNature nature_;  ///< Nature of the discrete field.

  Real west_bndy_;   ///< West boundary spatial coordinate.
  Real east_bndy_;   ///< East boundary spatial coordinate.
  int num_cells_x_;  ///< Number of cells discretizing the domain.
  Real delta_x_;     ///< Computed \f$ \Delta x \f$.

  Real south_bndy_;  ///< West boundary spatial coordinate.
  Real north_bndy_;  ///< East boundary spatial coordinate.
  int num_cells_y_;  ///< Number of cells discretizing the domain.
  Real delta_y_;     ///< Computed \f$ \Delta y \f$.
};
}
#endif  // End of: MTK_INCLUDE_UNI_STG_GRID_2D_H_
