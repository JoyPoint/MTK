/*!
\file mtk_uni_stg_grid_1d.h

\brief Definition of an 1D uniform staggered grid.

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

#ifndef MTK_INCLUDE_UNI_STG_GRID_1D_H_
#define MTK_INCLUDE_UNI_STG_GRID_1D_H_

#include <vector>
#include <string>

#include "mtk_roots.h"

namespace mtk {

/*!
\class UniStgGrid1D

\ingroup c06-grids

\brief Uniform 1D Staggered Grid.

Uniform 1D Staggered Grid.
*/
class UniStgGrid1D {
 public:
  /// \brief Prints the grid as a tuple of arrays.
  friend std::ostream& operator <<(std::ostream& stream, UniStgGrid1D &in);

  /// \brief Default constructor.
  UniStgGrid1D();

  /*!
  \brief Copy constructor.

  \param [in] grid Given grid.
  */
  UniStgGrid1D(const UniStgGrid1D &grid);

  /*!
  \brief Construct a grid based on spatial discretization parameters.

  \param[in] west_bndy_x Coordinate for the west boundary.
  \param[in] east_bndy_x Coordinate for the east boundary.
  \param[in] num_cells_x Number of cells of the required grid.
  \param[in] nature Nature of the discrete field to hold.

  \sa mtk::FieldNature
  */
  UniStgGrid1D(const Real &west_bndy_x,
               const Real &east_bndy_x,
               const int &num_cells_x,
               const mtk::FieldNature &nature = mtk::SCALAR);

  /// \brief Destructor.
  ~UniStgGrid1D();

  /*!
  \brief Provides access to west boundary spatial coordinate.

  \return West boundary spatial coordinate.
  */
  Real west_bndy_x() const;

  /*!
  \brief Provides access to east boundary spatial coordinate.

  \return East boundary spatial coordinate.
  */
  Real east_bndy_x() const;

  /*!
  \brief Provides access to the computed \$ \Delta x \$.

  \return Computed \$ \Delta x \$.
  */
  Real delta_x() const;

  /*!
  \brief Provides access to the grid spatial data.

  \return Pointer to the spatial data.

  \todo Review const-correctness of the pointer we return.
  */
  const Real *discrete_domain_x() const;

  /*!
  \brief Provides access to the grid field data.

  \return Pointer to the field data.

  \todo Review const-correctness of the pointer we return. Look at the STL!
  */
  Real *discrete_field();

  /*!
  \brief Provides access to the number of cells of the grid.

  \return Number of cells of the grid.
  */
  int num_cells_x() const;

  /*!
  \brief Binds a given scalar field to the grid.

  \param[in] ScalarField Pointer to the function implementing the scalar field.
  */
  void BindScalarField(Real (*ScalarField)(const Real &xx));

  /*!
  \brief Binds a given vector field to the grid.

  We assume the field to be of the form:
  \f[
    \mathbf{v}(\mathbf{x}) = v(x)\hat{\mathbf{i}}
  \f]

  \param[in] VectorField Pointer to the function implementing the vector field.
  */
  void BindVectorField(Real (*VectorField)(Real xx));

  /*!
  \brief Writes grid to a file compatible with gnuplot 4.6.

  \param[in] filename Name of the output file.
  \param[in] space_name Name for the first column of the data.
  \param[in] field_name Name for the second column of the data.

  \return Success of the file writing process.

  \sa http://www.gnuplot.info/
  */
  bool WriteToFile(std::string filename,
                   std::string space_name,
                   std::string field_name) const;

 private:
  FieldNature nature_;  ///< Nature of the discrete field.

  std::vector<Real> discrete_domain_x_; ///< Array of spatial data.
  std::vector<Real> discrete_field_;  ///< Array of field's data.

  Real west_bndy_x_;  ///< West boundary spatial coordinate.
  Real east_bndy_x_;  ///< East boundary spatial coordinate.
  Real num_cells_x_;  ///< Number of cells discretizing the domain.
  Real delta_x_;      ///< Produced \f$ \Delta x \f$.
};
}
#endif  // End of: MTK_INCLUDE_UNI_STG_GRID_1D_H_
