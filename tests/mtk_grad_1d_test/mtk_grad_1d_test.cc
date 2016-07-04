/*!
\file mtk_grad_1d_test.cc

\brief Testing the mimetic 1D gradient, constructed with the CBS algorithm.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
*/
/*
Copyright (C) 2016, Computational Science Research Center, San Diego State
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

#include <iostream>

#include "mtk.h"

void TestDefaultConstructorFactoryMethodDefault() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::Grad1D grad2;

  bool assertion = grad2.ConstructGrad1D();

  if (!assertion) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;

  }

  std::cout << grad2 << std::endl;

  mtk::Tools::EndUnitTestNo(1);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodFourthOrder() {

  mtk::Tools::BeginUnitTestNo(2);

  mtk::Grad1D grad4;

  bool assertion = grad4.ConstructGrad1D(4);

  if (!assertion) {
    std::cerr << "Mimetic grad (4th order) could not be built." << std::endl;
  }

  std::cout << grad4 << std::endl;

  mtk::Tools::EndUnitTestNo(2);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodSixthOrder() {

  mtk::Tools::BeginUnitTestNo(3);

  mtk::Grad1D grad6;

  bool assertion = grad6.ConstructGrad1D(6);

  if (!assertion) {
    std::cerr << "Mimetic grad (6th order) could not be built." << std::endl;
  }

  std::cout << grad6 << std::endl;

  mtk::Tools::EndUnitTestNo(3);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodEightOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(4);

  mtk::Grad1D grad8;

  bool assertion = grad8.ConstructGrad1D(8);

  if (!assertion) {
    std::cerr << "Mimetic grad (8th order) could not be built." << std::endl;
  }

  std::cout << grad8 << std::endl;

  mtk::Tools::EndUnitTestNo(4);
  mtk::Tools::Assert(assertion);
}

void TestDefaultConstructorFactoryMethodTenthOrderDefThreshold() {

  mtk::Tools::BeginUnitTestNo(5);

  mtk::Grad1D grad10;

  bool assertion = grad10.ConstructGrad1D(10);

  if (!assertion) {
    std::cerr << "Mimetic grad (10th order) could not be built." << std::endl;
  }

  std::cout << grad10 << std::endl;

  mtk::Tools::EndUnitTestNo(5);
  mtk::Tools::Assert(assertion);
}

void TestReturnAsDenseMatrixWithGrid() {

  mtk::Tools::BeginUnitTestNo(6);

  mtk::Grad1D grad2;

  bool assertion = grad2.ConstructGrad1D();

  if (!assertion) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;
  }

  mtk::UniStgGrid1D grid(0.0, 1.0, 5);

  mtk::DenseMatrix grad2m(grad2.ReturnAsDenseMatrix(grid));

  int rr{6};
  int cc{7};

  mtk::DenseMatrix ref(rr, cc);

  ref.set_encoded_operator(mtk::EncodedOperator::GRADIENT);

  // Row 1.
  ref.SetValue(0,0,-13.3333);
  ref.SetValue(0,1,15);
  ref.SetValue(0,2,-1.66667);
  ref.SetValue(0,3,0.0);
  ref.SetValue(0,4,0.0);
  ref.SetValue(0,5,0.0);
  ref.SetValue(0,6,0.0);

  // Row 2.
  ref.SetValue(1,0,0.0);
  ref.SetValue(1,1,-5.0);
  ref.SetValue(1,2,5.0);
  ref.SetValue(1,3,0.0);
  ref.SetValue(1,4,0.0);
  ref.SetValue(1,5,0.0);
  ref.SetValue(1,6,0.0);

  // Row 3.
  ref.SetValue(2,0,0.0);
  ref.SetValue(2,1,0.0);
  ref.SetValue(2,2,-5.0);
  ref.SetValue(2,3,5.0);
  ref.SetValue(2,4,0.0);
  ref.SetValue(2,5,0.0);
  ref.SetValue(2,6,0.0);

  // Row 4.
  ref.SetValue(3,0,0.0);
  ref.SetValue(3,1,0.0);
  ref.SetValue(3,2,0.0);
  ref.SetValue(3,3,-5.0);
  ref.SetValue(3,4,5.0);
  ref.SetValue(3,5,0.0);
  ref.SetValue(3,6,0.0);

  // Row 5.
  ref.SetValue(4,0,0.0);
  ref.SetValue(4,1,0.0);
  ref.SetValue(4,2,0.0);
  ref.SetValue(4,3,0.0);
  ref.SetValue(4,4,-5.0);
  ref.SetValue(4,5,5.0);
  ref.SetValue(4,6,0.0);

  // Row 6.
  ref.SetValue(5,0,0.0);
  ref.SetValue(5,1,0.0);
  ref.SetValue(5,2,0.0);
  ref.SetValue(5,3,0.0);
  ref.SetValue(5,4,1.66667);
  ref.SetValue(5,5,-15.0);
  ref.SetValue(5,6,13.3333);

  std::cout << grad2m << std::endl;
  std::cout << ref << std::endl;

  mtk::Tools::EndUnitTestNo(6);
  mtk::Tools::Assert(grad2m == ref);
}

void TestReturnAsDimensionlessDenseMatrix() {

  mtk::Tools::BeginUnitTestNo(7);

  mtk::Grad1D grad4;

  bool assertion = grad4.ConstructGrad1D(4);

  if (!assertion) {
    std::cerr << "Mimetic grad (4th order) could not be built." << std::endl;
  }

  mtk::DenseMatrix grad4m(grad4.ReturnAsDimensionlessDenseMatrix(10));

  std::cout << grad4m << std::endl;

  mtk::Tools::EndUnitTestNo(7);
  mtk::Tools::Assert(assertion);
}

void TestWriteToFile() {

  mtk::Tools::BeginUnitTestNo(8);

  mtk::Grad1D grad2;

  bool assertion = grad2.ConstructGrad1D();

  if (!assertion) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;
  }

  mtk::UniStgGrid1D grid(0.0, 1.0, 50);

  mtk::DenseMatrix grad2m(grad2.ReturnAsDenseMatrix(grid));

  std::cout << grad2m << std::endl;

  assertion = assertion && grad2m.WriteToFile("mtk_grad_1d_test_08_grad2m.dat");

  if(!assertion) {
    std::cerr << "Error writing to file." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(8);
  mtk::Tools::Assert(assertion);
}

void TestMimBndy() {

  mtk::Tools::BeginUnitTestNo(9);

  mtk::Grad1D grad2;

  bool assertion = grad2.ConstructGrad1D();

  if (!assertion) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;
  }

  std::cout << grad2 << std::endl;

  mtk::DenseMatrix grad2m(grad2.mim_bndy());

  std::cout << grad2m << std::endl;

  mtk::Tools::EndUnitTestNo(9);
  mtk::Tools::Assert(assertion);
}

int main () {

  std::cout << "Testing mtk::Grad1D class." << std::endl;

  TestDefaultConstructorFactoryMethodDefault();
  TestDefaultConstructorFactoryMethodFourthOrder();
  TestDefaultConstructorFactoryMethodSixthOrder();
  TestDefaultConstructorFactoryMethodEightOrderDefThreshold();
  TestDefaultConstructorFactoryMethodTenthOrderDefThreshold();
  TestReturnAsDenseMatrixWithGrid();
  TestReturnAsDimensionlessDenseMatrix();
  TestWriteToFile();
  TestMimBndy();
}
