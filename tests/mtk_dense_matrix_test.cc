/*!
\file mtk_dense_matrix_test.cc

\brief Test file for the mtk::DenseMatrix class.

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
#include <ctime>

#include "mtk.h"

void Test1() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::DenseMatrix m1;

  std::cout << m1 << std::endl;

  mtk::Tools::EndUnitTestNo(1);
}

void Test2() {

  mtk::Tools::BeginUnitTestNo(2);

  int rr = 4;
  int cc = 7;

  mtk::DenseMatrix m2(rr,cc);

  std::cout << m2 << std::endl;

  mtk::Tools::EndUnitTestNo(2);
}

void Test3() {

  mtk::Tools::BeginUnitTestNo(3);

  int rank = 5;
  bool padded = true;
  bool transpose = false;

  mtk::DenseMatrix m3(rank,padded,transpose);

  std::cout << m3 << std::endl;

  mtk::Tools::EndUnitTestNo(3);
}

void Test4() {

  mtk::Tools::BeginUnitTestNo(4);

  int rank = 5;
  bool padded = false;
  bool transpose = false;

  mtk::DenseMatrix m4(rank,padded,transpose);

  std::cout << m4 << std::endl;

  mtk::Tools::EndUnitTestNo(4);
}

void Test5() {

  mtk::Tools::BeginUnitTestNo(5);

  int rr = 4;
  int cc = 7;

  mtk::DenseMatrix m5(rr,cc);

  for (auto ii = 0; ii < rr; ++ii) {
    for (auto jj = 0; jj < cc; ++jj) {
      m5.SetValue(ii,jj,(mtk::Real) ii + jj);
    }
  }

  std::cout << m5 << std::endl;

  mtk::Real *vals = m5.data();

  for (auto ii = 0; ii < rr; ++ii) {
    for (auto jj = 0; jj < cc; ++jj) {
      std::cout << " " << vals[ii*cc + jj];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for (auto ii = 0; ii < rr; ++ii) {
    for (auto jj = 0; jj < cc; ++jj) {
      std::cout << " " << m5.GetValue(ii,jj);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  mtk::Tools::EndUnitTestNo(5);
}

void Test6() {

  mtk::Tools::BeginUnitTestNo(6);

  bool transpose = false;
  int generator_length = 3;
  int progression_length = 4;

  mtk::Real generator[] = {-0.5, 0.5, 1.5};

  mtk::DenseMatrix m6(generator,generator_length,progression_length,transpose);

  std::cout << m6 << std::endl;

  transpose = true;

  mtk::DenseMatrix m7(generator,generator_length,progression_length,transpose);

  std::cout << m7 << std::endl;


  mtk::Tools::EndUnitTestNo(6);
}

void Test7() {

  mtk::Tools::BeginUnitTestNo(7);

  bool padded = false;
  bool transpose = false;
  int lots_of_rows = 2;
  int lots_of_cols = 5;
  mtk::DenseMatrix m8(lots_of_rows,padded,transpose);

  std::cout << m8 << std::endl;

  mtk::DenseMatrix m9(lots_of_rows,lots_of_cols);

  for (auto ii = 0; ii < lots_of_rows; ++ii) {
    for (auto jj = 0; jj < lots_of_cols; ++jj) {
      m9.SetValue(ii,jj,(mtk::Real) ii*lots_of_cols + jj + 1);
    }
  }

  std::cout << m9 << std::endl;

  mtk::DenseMatrix m10 = mtk::DenseMatrix::Kron(m8,m9);

  std::cout << m10 << std::endl;

  mtk::Tools::EndUnitTestNo(7);
}

void Test8() {

  mtk::Tools::BeginUnitTestNo(8);

  int lots_of_rows = 4;
  int lots_of_cols = 3;
  mtk::DenseMatrix m11(lots_of_rows,lots_of_cols);

  for (auto ii = 0; ii < lots_of_rows; ++ii) {
    for (auto jj = 0; jj < lots_of_cols; ++jj) {
      m11.SetValue(ii,jj,(mtk::Real) ii*lots_of_cols + jj + 1);
    }
  }

  std::cout << m11 << std::endl;

  m11.Transpose();

  std::cout << m11 << std::endl;

  mtk::DenseMatrix m12;

  m12 = m11;

  std::cout << m12 << std::endl;

  mtk::Tools::EndUnitTestNo(8);
}

void Test9() {

  mtk::Tools::BeginUnitTestNo(9);

  bool transpose = false;
  int gg_l = 3;
  int progression_length = 4;
  mtk::Real gg[] = {-0.5, 0.5, 1.5};

  mtk::DenseMatrix m13(gg, gg_l ,progression_length, transpose);

  std::cout << m13 << std::endl;

  mtk::DenseMatrix m14;

  m14 = m13;

  std::cout << m14 << std::endl;

  m13.Transpose();

  std::cout << m13 << std::endl;

  m14 = m13;

  std::cout << m14 << std::endl;

  mtk::Tools::EndUnitTestNo(9);
}

int main () {

  std::cout << "Testing mtk::DenseMatrix class." << std::endl;

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
