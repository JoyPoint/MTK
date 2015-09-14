#include <vector>

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"
#include "mtk_uni_stg_grid_1d.h"

namespace mtk {

class BCDesc1D {
 public:
  static void ImposeOnOperator(DenseMatrix &matrix,
                               const std::vector<Real> &west,
                               const std::vector<Real> &east);

  static void ImposeOnGrid(UniStgGrid1D &grid,
                           const Real &omega,
                           const Real &epsilon);
};
}
