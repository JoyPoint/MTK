#include "mtk_tools.h"

#include "mtk_bc_desc_1d.h"

void mtk::BCDesc1D::ImposeOnOperator(mtk::DenseMatrix &matrix,
                                     const std::vector<mtk::Real> &west,
                                     const std::vector<mtk::Real> &east) {

  mtk::Tools::Prevent(matrix.num_rows() == 0, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(west.size() > (unsigned int) matrix.num_cols(),
                      __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(east.size() > (unsigned int) matrix.num_cols(),
                      __FILE__, __LINE__, __func__);

  /// 1. Assign the west array.

  for (unsigned int ii = 0; ii < west.size(); ++ii) {
    matrix.SetValue(0, ii, west[ii]);
  }

  /// 2. Assign the east array.

  for (unsigned int ii = 0; ii < east.size(); ++ii) {
    matrix.SetValue(matrix.num_rows() - 1,
                    matrix.num_cols() - 1 - ii,
                    east[ii]);
  }
}

void mtk::BCDesc1D::ImposeOnGrid(mtk::UniStgGrid1D &grid,
                                 const mtk::Real &omega,
                                 const mtk::Real &epsilon) {

  mtk::Tools::Prevent(grid.num_cells_x() == 0, __FILE__, __LINE__, __func__);

  /// 1. Assign the west condition.

  grid.discrete_field_u()[0] = omega;

  /// 2. Assign the east condition.

  grid.discrete_field_u()[grid.num_cells_x() + 2 - 1] = epsilon;
}
