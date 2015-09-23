#ifndef MTK_INCLUDE_MTK_GRAD_2D_H_
#define MTK_INCLUDE_MTK_GRAD_2D_H_

namespace mtk{

class Grad2D {
 public:
  friend std::ostream &operator<<(std::ostream& stream, const Grad2D &grad);
  Grad2D();
  Grad2D(const Grad2D &grad);
  ~Grad2D();
  bool ConstructGrad2D(int order_accuracy = kDefaultOrderAccuracy,
                       Real mimetic_threshold = kDefaultMimeticThreshold);
  DenseMatrix ReturnAsDenseMatrix(const UniStgGrid2D &grid);

 private:
  int order_accuracy_;
  Real mimetic_threshold_;
};
}
#endif  // End of: MTK_INCLUDE_MTK_GRAD_2D_H_
