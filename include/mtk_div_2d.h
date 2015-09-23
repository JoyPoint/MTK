#ifndef MTK_INCLUDE_MTK_DIV_2D_H_
#define MTK_INCLUDE_MTK_DIV_2D_H_

namespace mtk{

class Div2D {
 public:
  friend std::ostream &operator<<(std::ostream& stream, const Div2D &div);
  Div2D();
  Div2D(const Div2D &div);
  ~Div2D();
  bool ConstructDiv2D(int order_accuracy = kDefaultOrderAccuracy,
                      Real mimetic_threshold = kDefaultMimeticThreshold);
  DenseMatrix ReturnAsDenseMatrix(const UniStgGrid2D &grid);

 private:
  int order_accuracy_;
  Real mimetic_threshold_;
};
}
#endif  // End of: MTK_INCLUDE_MTK_DIV_2D_H_
