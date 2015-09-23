#ifndef MTK_INCLUDE_MTK_INTERP_2D_H_
#define MTK_INCLUDE_MTK_INTERP_2D_H_

namespace mtk{

class Interp2D {
 public:
  friend std::ostream  &operator<<(std::ostream& stream, const Interp2D &inter);
  Interp2D();
  Interp2D(const Interp2D &inter);
  ~Interp2D();
  bool ConstructInterp2D(int order_accuracy = kDefaultOrderAccuracy,
                         Real mimetic_threshold = kDefaultMimeticThreshold);
  DenseMatrix ReturnAsDenseMatrix(const UniStgGrid2D &grid);

 private:
  int order_accuracy_;
  Real mimetic_threshold_;
};
}
#endif  // End of: MTK_INCLUDE_MTK_INTERP_2D_H_
