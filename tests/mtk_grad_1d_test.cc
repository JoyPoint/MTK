/*!
\file mtk_grad_1d_test.cc

\brief Testing the mimetic 1D gradient, constructed with the CBS algorithm.

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

  mtk::Tools::BeginTestNo(1);

  mtk::Grad1D grad2;

  bool info = grad2.ConstructGrad1D();

  if (!info) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;
  }

  std::cout << grad2 << std::endl;

  mtk::Tools::EndTestNo(1);
}

void Test2() {

  mtk::Tools::BeginTestNo(2);

  mtk::Grad1D grad4;

  bool info = grad4.ConstructGrad1D(4);

  if (!info) {
    std::cerr << "Mimetic grad (4th order) could not be built." << std::endl;
  }

  std::cout << grad4 << std::endl;

  mtk::Tools::EndTestNo(2);
}

void Test3() {

  mtk::Tools::BeginTestNo(3);

  mtk::Grad1D grad6;

  bool info = grad6.ConstructGrad1D(6);

  if (!info) {
    std::cerr << "Mimetic grad (6th order) could not be built." << std::endl;
  }

  std::cout << grad6 << std::endl;

  mtk::Tools::EndTestNo(3);
}

void Test4() {

  mtk::Tools::BeginTestNo(4);

  mtk::Grad1D grad8;

  bool info = grad8.ConstructGrad1D(8);

  if (!info) {
    std::cerr << "Mimetic grad (8th order) could not be built." << std::endl;
  }

  std::cout << grad8 << std::endl;

  mtk::Tools::EndTestNo(4);
}

void Test5() {

  mtk::Tools::BeginTestNo(5);

  mtk::Grad1D grad10;

  bool info = grad10.ConstructGrad1D(10);

  if (!info) {
    std::cerr << "Mimetic grad (10th order) could not be built." << std::endl;
  }

  std::cout << grad10 << std::endl;

  mtk::Tools::EndTestNo(5);
}

void Test6() {

  mtk::Tools::BeginTestNo(6);

  mtk::Grad1D grad2;

  bool info = grad2.ConstructGrad1D();

  if (!info) {
    std::cerr << "Mimetic grad (2nd order) could not be built." << std::endl;
  }

  std::cout << grad2 << std::endl;

  mtk::UniStgGrid1D grid(0.0, 1.0, 5);

  std::cout << grid << std::endl;

  mtk::DenseMatrix grad2m(grad2.ReturnAsDenseMatrix(grid));

  std::cout << grad2m << std::endl;

  mtk::Tools::EndTestNo(6);
}

int main () {

  std::cout << "Testing mtk::Grad1D class." << std::endl;

  Test1();
  Test2();
  Test3();
  Test4();
  Test5();
  Test6();
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
