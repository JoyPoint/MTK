#ifndef MTK_INCLUDE_UNI_STG_GRID_2D_H_
#define MTK_INCLUDE_UNI_STG_GRID_2D_H_

#include <vector>
#include <string>

#include "mtk_roots.h"
#include "mtk_enums.h"

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
  \brief Provides access to the computed \$ \Delta x \$.

  \return Computed \$ \Delta x \$.
  */
  Real delta_x() const;

  /*!
  \brief Provides access to the computed \$ \Delta y \$.

  \return Computed \$ \Delta y \$.
  */
  Real delta_y() const;

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
  Real *discrete_field_u();

  /*!
  \brief Provides access to the number of cells of the grid.

  \return Number of cells of the grid.
  */
  int num_cells_x() const;

  /*!
  \brief Provides access to the number of cells of the grid.

  \return Number of cells of the grid.
  */
  int num_cells_y() const;

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

  \param[in] VectorField Pointer to the function implementing the vector field.
  */
  void BindVectorFieldPComponent(Real (*VectorField)(Real xx, Real yy));

  /*!
  \brief Binds a given vector field to the grid.

  We assume the field to be of the form:

  \f[
    \mathbf{v}(x) = p(x, y)\hat{\mathbf{i}} + q(x, y)\hat{\mathbf{j}}
  \f]

  \param[in] VectorField Pointer to the function implementing the vector field.
  */
  void BindVectorFieldQComponent(Real (*VectorField)(Real xx, Real yy));

  /*!
  \brief Writes grid to a file compatible with Gnuplot 4.6.

  \param[in] filename Name of the output file.
  \param[in] space_name Name for the first column of the data.
  \param[in] field_name Name for the second column of the data.

  \return Success of the file writing process.

  \sa http://www.gnuplot.info/
  */
  bool WriteToFile(std::string filename,
                   std::string space_name,
                   std::string field_name);

 private:
  FieldNature nature_;  ///< Nature of the discrete field.

  std::vector<Real> discrete_domain_x_; ///< Array of spatial data.
  std::vector<Real> discrete_domain_y_; ///< Array of spatial data.
  std::vector<Real> discrete_field_u_;  ///< Array of field's data.

  Real west_bndy_x_;  ///< West boundary spatial coordinate.
  Real east_bndy_x_;  ///< East boundary spatial coordinate.
  Real num_cells_x_;  ///< Number of cells discretizing the domain.
  Real delta_x_;      ///< Produced \f$ \Delta x\f$.

  Real south_bndy_y_;  ///< West boundary spatial coordinate.
  Real north_bndy_y_;  ///< East boundary spatial coordinate.
  Real num_cells_y_;  ///< Number of cells discretizing the domain.
  Real delta_y_;      ///< Produced \f$ \Delta y\f$.
};
}
#endif  // End of: MTK_INCLUDE_UNI_STG_GRID_2D_H_
