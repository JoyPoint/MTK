/*!
\file mtk_div_1d_test.cc

\brief Testing the mimetic 1D divergence, constructed with the CBS algorithm.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
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

#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"

void TestDefaultConstructorFactoryMethodDefault() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::Div1D div2;

  bool assertion = div2.ConstructDiv1D();

  if (!assertion) {
    std::cerr << "Mimetic div (2nd order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(1);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodFourthOrder() {

  mtk::Tools::BeginUnitTestNo(2);

  mtk::Div1D div4;

  bool assertion = div4.ConstructDiv1D(4);

  if (!assertion) {
    std::cerr << "Mimetic div (4th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(2);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodSixthOrder() {

  mtk::Tools::BeginUnitTestNo(3);

  mtk::Div1D div6;

  bool assertion = div6.ConstructDiv1D(6);

  if (!assertion) {
    std::cerr << "Mimetic div (6th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(3);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodEightOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(4);

  mtk::Div1D div8;

  bool assertion = div8.ConstructDiv1D(8);

  if (!assertion) {
    std::cerr << "Mimetic div (8th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(4);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodTenthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(5);

  mtk::Div1D div10;

  bool assertion = div10.ConstructDiv1D(10);

  if (!assertion) {
    std::cerr << "Mimetic div (10th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(5);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodTwelfthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(6);

  mtk::Div1D div12;

  bool assertion = div12.ConstructDiv1D(12);

  if (!assertion) {
    std::cerr << "Mimetic div (12th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(6);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodFourteenthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(7);

  mtk::Div1D div14;

  bool assertion = div14.ConstructDiv1D(14);

  if (!assertion) {
    std::cerr << "Mimetic div (14th order) could not be built." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(7);
  mtk::Tools::Assert(assertion);
}

void TestSecondOrderReturnAsDenseMatrixWithGrid() {

  mtk::Tools::BeginUnitTestNo(8);

  mtk::Div1D div2;

  bool assertion = div2.ConstructDiv1D();

  if (!assertion) {
    std::cerr << "Mimetic div (2nd order) could not be built." << std::endl;
  }

  mtk::UniStgGrid1D grid(0.0, 1.0, 5);

  mtk::DenseMatrix div2m(div2.ReturnAsDenseMatrix(grid));

  int rr{7};
  int cc{6};

  mtk::DenseMatrix ref(rr, cc);

  // Row 2.
  ref.SetValue(1,0,-5.0);
  ref.SetValue(1,1,5.0);
  ref.SetValue(1,2,0.0);
  ref.SetValue(1,3,0.0);
  ref.SetValue(1,4,0.0);
  ref.SetValue(1,5,0.0);
  ref.SetValue(1,6,0.0);

  // Row 3.
  ref.SetValue(2,0,0.0);
  ref.SetValue(2,1,-5.0);
  ref.SetValue(2,2,5.0);
  ref.SetValue(2,3,0.0);
  ref.SetValue(2,4,0.0);
  ref.SetValue(2,5,0.0);
  ref.SetValue(2,6,0.0);

  // Row 4.
  ref.SetValue(3,0,0.0);
  ref.SetValue(3,1,0.0);
  ref.SetValue(3,2,-5.0);
  ref.SetValue(3,3,5.0);
  ref.SetValue(3,4,0.0);
  ref.SetValue(3,5,0.0);
  ref.SetValue(3,6,0.0);

  // Row 5.
  ref.SetValue(4,0,0.0);
  ref.SetValue(4,1,0.0);
  ref.SetValue(4,2,0.0);
  ref.SetValue(4,3,-5.0);
  ref.SetValue(4,4,5.0);
  ref.SetValue(4,5,0.0);
  ref.SetValue(4,6,0.0);

  // Row 6.
  ref.SetValue(5,0,0.0);
  ref.SetValue(5,1,0.0);
  ref.SetValue(5,2,0.0);
  ref.SetValue(5,3,0.0);
  ref.SetValue(5,4,-5.0);
  ref.SetValue(5,5,5.0);
  ref.SetValue(5,6,0.0);

  assertion = assertion && (div2m == ref);

  mtk::Tools::EndUnitTestNo(8);
  mtk::Tools::Assert(assertion);
}

void TestFourthOrderReturnAsDenseMatrixWithGrid() {

  mtk::Tools::BeginUnitTestNo(9);

  mtk::Div1D div4;

  bool assertion = div4.ConstructDiv1D(4);

  if (!assertion) {
    std::cerr << "Mimetic div (4th order) could not be built." << std::endl;
  }

  std::cout << div4 << std::endl;

  mtk::UniStgGrid1D grid(0.0, 1.0, 11);

  std::cout << grid << std::endl;

  mtk::DenseMatrix div4m(div4.ReturnAsDenseMatrix(grid));

  std::cout << div4m << std::endl;

  mtk::Tools::EndUnitTestNo(9);
}

int main () {

  std::cout << "Testing mtk::Div1D class." << std::endl;

  TestDefaultConstructorFactoryMethodDefault();
  TestDefaultConstructorFactoryMethodFourthOrder();
  TestDefaultConstructorFactoryMethodSixthOrder();
  TestDefaultConstructorFactoryMethodEightOrderDefThreshold();
  TestDefaultConstructorFactoryMethodTenthOrderDefThreshold();
  TestDefaultConstructorFactoryMethodTwelfthOrderDefThreshold();
  TestDefaultConstructorFactoryMethodFourteenthOrderDefThreshold();
  TestSecondOrderReturnAsDenseMatrixWithGrid();
  TestFourthOrderReturnAsDenseMatrixWithGrid();
}

#else
#include <iostream>
using std::cout;
using std::endl;
int main () {
  cout << "This code HAS to be compiled with support for C++11." << endl;
  cout << "Exiting..." << endl;
}
#endif
