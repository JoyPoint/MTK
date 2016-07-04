/*!
\file mtk_lap_1d_test.cc

\brief Testing the 1D Laplacian operator.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\author: Johnny Corbino - jcorbino at mail dot sdsu dot edu
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
