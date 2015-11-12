#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"

void TestDefaultConstructorFactoryMethodDefault() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::Lap1D lap2;

  bool assertion = lap2.ConstructLap1D();

  if (!assertion) {
    std::cerr << "Mimetic lap (2nd order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(1);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodFourthOrder() {

  mtk::Tools::BeginUnitTestNo(2);

  mtk::Lap1D lap4;

  bool assertion = lap4.ConstructLap1D(4);

  if (!assertion) {
    std::cerr << "Mimetic lap (4th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(2);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodSixthOrder() {

  mtk::Tools::BeginUnitTestNo(3);

  mtk::Lap1D lap6;

  bool assertion = lap6.ConstructLap1D(6);

  if (!assertion) {
    std::cerr << "Mimetic lap (6th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(3);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodEightOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(4);

  mtk::Lap1D lap8;

  bool assertion = lap8.ConstructLap1D(8);

  if (!assertion) {
    std::cerr << "Mimetic lap (8th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(4);
}

void TestDefaultConstructorFactoryMethodTenthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(5);

  mtk::Lap1D lap10;

  bool assertion = lap10.ConstructLap1D(10);

  if (!assertion) {
    std::cerr << "Mimetic lap (10th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(5);
}

void TestDefaultConstructorFactoryMethodTwelfthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(6);

  mtk::Lap1D lap12;

  bool assertion = lap12.ConstructLap1D(12);

  if (!assertion) {
    std::cerr << "Mimetic lap (12th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(6);
}

void TestReturnAsDenseMatrix() {

  mtk::Tools::BeginUnitTestNo(8);

  mtk::Lap1D lap4;

  bool assertion = lap4.ConstructLap1D(4);

  if (!assertion) {
    std::cerr << "Mimetic lap (4th order) could not be built." << std::endl;
  }

  mtk::UniStgGrid1D aux(0.0, 1.0, 11);

  mtk::DenseMatrix lap4_m(lap4.ReturnAsDenseMatrix(aux));

  assertion = assertion &&
      abs(lap4_m.GetValue(1, 0) - 385.133) < mtk::kDefaultTolerance &&
      abs(lap4_m.GetValue(11, 12) - 385.133) < mtk::kDefaultTolerance;
  mtk::Tools::EndUnitTestNo(8);
  mtk::Tools::Assert(assertion);
}

int main () {

  std::cout << "Testing MTK 1D Laplacian" << std::endl;
  
  TestDefaultConstructorFactoryMethodDefault();
  TestDefaultConstructorFactoryMethodFourthOrder();
  TestDefaultConstructorFactoryMethodSixthOrder();
  TestDefaultConstructorFactoryMethodEightOrderDefThreshold();
  TestDefaultConstructorFactoryMethodTenthOrderDefThreshold();
  TestDefaultConstructorFactoryMethodTwelfthOrderDefThreshold();
  TestReturnAsDenseMatrix();
}

#else
#include <iostream>
int main () {
  std::cout << "This code HAS to be compiled to support C++11." << std::endl;
  std::cout << "Exiting..." << std::endl;
}
#endif
