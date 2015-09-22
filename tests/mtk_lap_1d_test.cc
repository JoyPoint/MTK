#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"

void Test1() {

  mtk::Tools::BeginTestNo(1);

  mtk::Lap1D lap2;

  bool info = lap2.ConstructLap1D();

  if (!info) {
    std::cerr << "Mimetic lap (2nd order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(1);
}

void Test2() {

  mtk::Tools::BeginTestNo(2);

  mtk::Lap1D lap4;

  bool info = lap4.ConstructLap1D(4);

  if (!info) {
    std::cerr << "Mimetic lap (4th order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(2);
}

void Test3() {

  mtk::Tools::BeginTestNo(3);

  mtk::Lap1D lap6;

  bool info = lap6.ConstructLap1D(6);

  if (!info) {
    std::cerr << "Mimetic lap (6th order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(3);
}

void Test4() {

  mtk::Tools::BeginTestNo(4);

  mtk::Lap1D lap8;

  bool info = lap8.ConstructLap1D(8);

  if (!info) {
    std::cerr << "Mimetic lap (8th order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(4);
}

void Test5() {

  mtk::Tools::BeginTestNo(5);

  mtk::Lap1D lap10;

  bool info = lap10.ConstructLap1D(10);

  if (!info) {
    std::cerr << "Mimetic lap (10th order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(5);
}

void Test6() {

  mtk::Tools::BeginTestNo(6);

  mtk::Lap1D lap12;

  bool info = lap12.ConstructLap1D(12);

  if (!info) {
    std::cerr << "Mimetic lap (12th order) could not be built." << std::endl;
  }

  mtk::Tools::EndTestNo(6);
}

void Test7() {

  mtk::Tools::BeginTestNo(7);

  mtk::Lap1D lap4;

  bool info = lap4.ConstructLap1D(4);

  if (!info) {
    std::cerr << "Mimetic lap (4th order) could not be built." << std::endl;
  }

  std::cout << lap4 << std::endl;
  std::cout << std::endl;

  mtk::Tools::EndTestNo(7);
}

void Test8() {

  mtk::Tools::BeginTestNo(8);

  mtk::Lap1D lap4;

  bool info = lap4.ConstructLap1D(4);

  if (!info) {
    std::cerr << "Mimetic lap (4th order) could not be built." << std::endl;
  }

  std::cout << lap4 << std::endl;
  std::cout << std::endl;

  mtk::UniStgGrid1D aux(0.0, 1.0, 11);

  mtk::DenseMatrix lap4_m(lap4.ReturnAsDenseMatrix(aux));

  std::cout << lap4_m << std::endl;
  std::cout << std::endl;

  mtk::Tools::EndTestNo(8);
}

int main () {

  std::cout << "Testing MTK 1D Laplacian" << std::endl;
  
  Test1();
  Test2();
  Test3();
  Test4();
  Test5();
  Test6();
  Test7();
  Test8();
}

#else
#include <iostream>
int main () {
  std::cout << "This code HAS to be compiled to support C++11." << std::endl;
  std::cout << "Exiting..." << std::endl;
}
#endif
