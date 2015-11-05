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

void Test1() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::Div1D div2;

  bool info = div2.ConstructDiv1D();

  if (!info) {
    std::cerr << "Mimetic div (2nd order) could not be built." << std::endl;
  }

  std::cout << div2 << std::endl;

  mtk::Tools::EndUnitTestNo(1);
}

void Test2() {

  mtk::Tools::BeginUnitTestNo(2);

  mtk::Div1D div4;

  bool info = div4.ConstructDiv1D(4);

  if (!info) {
    std::cerr << "Mimetic div (4th order) could not be built." << std::endl;
  }

  std::cout << div4 << std::endl;

  mtk::Tools::EndUnitTestNo(2);
}

void Test3() {

  mtk::Tools::BeginUnitTestNo(3);

  mtk::Div1D div6;

  bool info = div6.ConstructDiv1D(6);

  if (!info) {
    std::cerr << "Mimetic div (6th order) could not be built." << std::endl;
  }

  std::cout << div6 << std::endl;

  mtk::Tools::EndUnitTestNo(3);
}

void Test4() {

  mtk::Tools::BeginUnitTestNo(4);

  mtk::Div1D div8;

  bool info = div8.ConstructDiv1D(8);

  if (!info) {
    std::cerr << "Mimetic div (8th order) could not be built." << std::endl;
  }

  std::cout << div8 << std::endl;

  mtk::Tools::EndUnitTestNo(4);
}

void Test5() {

  mtk::Tools::BeginUnitTestNo(5);

  mtk::Div1D div10;

  bool info = div10.ConstructDiv1D(10);

  if (!info) {
    std::cerr << "Mimetic div (10th order) could not be built." << std::endl;
  }

  std::cout << div10 << std::endl;

  mtk::Tools::EndUnitTestNo(5);
}

void Test6() {

  mtk::Tools::BeginUnitTestNo(6);

  mtk::Div1D div12;

  bool info = div12.ConstructDiv1D(12);

  if (!info) {
    std::cerr << "Mimetic div (12th order) could not be built." << std::endl;
  }

  std::cout << div12 << std::endl;

  mtk::Tools::EndUnitTestNo(6);
}

void Test7() {

  mtk::Tools::BeginUnitTestNo(7);

  mtk::Div1D div14;

  bool info = div14.ConstructDiv1D(14);

  if (!info) {
    std::cerr << "Mimetic div (14th order) could not be built." << std::endl;
  }

  std::cout << div14 << std::endl;

  mtk::Tools::EndUnitTestNo(7);
}

void Test8() {

  mtk::Tools::BeginUnitTestNo(8);

  mtk::Div1D div2;

  bool info = div2.ConstructDiv1D();

  if (!info) {
    std::cerr << "Mimetic div (2nd order) could not be built." << std::endl;
  }

  std::cout << div2 << std::endl;

  mtk::UniStgGrid1D grid(0.0, 1.0, 5);

  std::cout << grid << std::endl;

  mtk::DenseMatrix div2m(div2.ReturnAsDenseMatrix(grid));

  std::cout << div2m << std::endl;

  mtk::Tools::EndUnitTestNo(8);
}

void Test9() {

  mtk::Tools::BeginUnitTestNo(9);

  mtk::Div1D div4;

  bool info = div4.ConstructDiv1D(4);

  if (!info) {
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

  Test1();
  Test2();
  Test3();
  Test4();
  Test5();
  Test6();
  Test7();
  Test8();
  Test9();
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
