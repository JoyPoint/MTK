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

#if __cplusplus == 201103L

#include <iostream>
#include <ctime>

#include "mtk.h"

void TestDefaultConstructor() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::DenseMatrix m1;

  mtk::Tools::EndUnitTestNo(1);
  mtk::Tools::Assert(m1.data() == nullptr);
}

void TestConstructorWithNumRowsNumCols() {

  mtk::Tools::BeginUnitTestNo(2);

  int rr = 4;
  int cc = 7;

  mtk::DenseMatrix m2(rr,cc);

  mtk::Tools::EndUnitTestNo(2);

  bool assertion =
    m2.data() != nullptr && m2.num_rows() == rr && m2.num_cols() == cc;

  mtk::Tools::Assert(assertion);
}

void TestConstructAsIdentity() {

  mtk::Tools::BeginUnitTestNo(3);

  int rank = 5;
  bool padded = true;
  bool transpose = false;

  mtk::DenseMatrix m3(rank,padded,transpose);

  mtk::DenseMatrix rr(rank + 2,rank);

  for (int ii = 0; ii < rank; ++ii) {
    rr.SetValue(ii + 1, ii, mtk::kOne);
  }

  mtk::Tools::EndUnitTestNo(3);
  mtk::Tools::Assert(m3 == rr);
}

void TestConstructAsVandermonde() {

  mtk::Tools::BeginUnitTestNo(4);

  int rank = 5;
  bool padded = false;
  bool transpose = false;

  mtk::DenseMatrix m4(rank,padded,transpose);

  mtk::DenseMatrix rr(rank,rank);

  for (int ii = 0; ii < rank; ++ii) {
    rr.SetValue(ii, ii, mtk::kOne);
  }

  mtk::Tools::EndUnitTestNo(4);
  mtk::Tools::Assert(m4 == rr);
}

void TestSetValueGetValue() {

  mtk::Tools::BeginUnitTestNo(5);

  int rr = 4;
  int cc = 7;

  mtk::DenseMatrix m5(rr,cc);

  for (auto ii = 0; ii < rr; ++ii) {
    for (auto jj = 0; jj < cc; ++jj) {
      m5.SetValue(ii,jj,(mtk::Real) ii + jj);
    }
  }

  mtk::Real *vals = m5.data();

  bool assertion{true};

  for (auto ii = 0; ii < rr && assertion; ++ii) {
    for (auto jj = 0; jj < cc  && assertion; ++jj) {
      assertion = assertion && m5.GetValue(ii,jj) == vals[ii*cc + jj];
    }
  }

  mtk::Tools::EndUnitTestNo(5);
  mtk::Tools::Assert(assertion);
}

void TestConstructAsVandermondeTranspose() {

  mtk::Tools::BeginUnitTestNo(6);

  bool transpose = false;
  int generator_length = 3;
  int progression_length = 4;

  mtk::Real generator[] = {-0.5, 0.5, 1.5};

  mtk::DenseMatrix m6(generator,generator_length,progression_length,transpose);

  transpose = true;

  mtk::DenseMatrix m7(generator,generator_length,progression_length,transpose);
  mtk::DenseMatrix rr(progression_length, generator_length);

  rr.SetValue(0, 0, 1.0);
  rr.SetValue(0, 1, 1.0);
  rr.SetValue(0, 2, 1.0);

  rr.SetValue(1, 0, -0.5);
  rr.SetValue(1, 1, 0.5);
  rr.SetValue(1, 2, 1.5);

  rr.SetValue(2, 0, 0.25);
  rr.SetValue(2, 1, 0.25);
  rr.SetValue(2, 2, 2.25);

  rr.SetValue(3, 0, -0.125);
  rr.SetValue(3, 1, 0.125);
  rr.SetValue(3, 2, 3.375);

  mtk::Tools::EndUnitTestNo(6);
  mtk::Tools::Assert(m7 == rr);
}

void TestKron() {

  mtk::Tools::BeginUnitTestNo(7);

  bool padded = false;
  bool transpose = false;
  int lots_of_rows = 2;
  int lots_of_cols = 5;
  mtk::DenseMatrix m8(lots_of_rows,padded,transpose);

  mtk::DenseMatrix m9(lots_of_rows,lots_of_cols);

  for (auto ii = 0; ii < lots_of_rows; ++ii) {
    for (auto jj = 0; jj < lots_of_cols; ++jj) {
      m9.SetValue(ii,jj,(mtk::Real) ii*lots_of_cols + jj + 1);
    }
  }

  mtk::DenseMatrix m10 = mtk::DenseMatrix::Kron(m8,m9);

  mtk::DenseMatrix rr(lots_of_rows*lots_of_rows, lots_of_rows*lots_of_cols);

  rr.SetValue(0,0,1.0);
  rr.SetValue(0,1,2.0);
  rr.SetValue(0,2,3.0);
  rr.SetValue(0,3,4.0);
  rr.SetValue(0,4,5.0);
  rr.SetValue(0,5,0.0);
  rr.SetValue(0,6,0.0);
  rr.SetValue(0,7,0.0);
  rr.SetValue(0,8,0.0);
  rr.SetValue(0,9,0.0);

  rr.SetValue(1,0,6.0);
  rr.SetValue(1,1,7.0);
  rr.SetValue(1,2,8.0);
  rr.SetValue(1,3,9.0);
  rr.SetValue(1,4,10.0);
  rr.SetValue(1,5,0.0);
  rr.SetValue(1,6,0.0);
  rr.SetValue(1,7,0.0);
  rr.SetValue(1,8,0.0);
  rr.SetValue(1,9,0.0);

  rr.SetValue(2,0,0.0);
  rr.SetValue(2,1,0.0);
  rr.SetValue(2,2,0.0);
  rr.SetValue(2,3,0.0);
  rr.SetValue(2,4,0.0);
  rr.SetValue(2,5,1.0);
  rr.SetValue(2,6,2.0);
  rr.SetValue(2,7,3.0);
  rr.SetValue(2,8,4.0);
  rr.SetValue(2,9,5.0);

  rr.SetValue(3,0,0.0);
  rr.SetValue(3,1,0.0);
  rr.SetValue(3,2,0.0);
  rr.SetValue(3,3,0.0);
  rr.SetValue(3,4,0.0);
  rr.SetValue(3,5,6.0);
  rr.SetValue(3,6,7.0);
  rr.SetValue(3,7,8.0);
  rr.SetValue(3,8,9.0);
  rr.SetValue(3,9,10.0);

  mtk::Tools::EndUnitTestNo(7);
  mtk::Tools::Assert(m10 == rr);
}

void TestConstructWithNumRowsNumColsAssignmentOperator() {

  mtk::Tools::BeginUnitTestNo(8);

  int lots_of_rows = 4;
  int lots_of_cols = 3;
  mtk::DenseMatrix m11(lots_of_rows,lots_of_cols);

  for (auto ii = 0; ii < lots_of_rows; ++ii) {
    for (auto jj = 0; jj < lots_of_cols; ++jj) {
      m11.SetValue(ii,jj,(mtk::Real) ii*lots_of_cols + jj + 1);
    }
  }

  m11.Transpose();

  mtk::DenseMatrix m12;

  m12 = m11;

  mtk::Tools::EndUnitTestNo(8);
  mtk::Tools::Assert(m11 == m12);
}

void TestConstructAsVandermondeTransposeAssignmentOperator() {

  mtk::Tools::BeginUnitTestNo(9);

  bool transpose = false;
  int gg_l = 3;
  int progression_length = 4;
  mtk::Real gg[] = {-0.5, 0.5, 1.5};

  mtk::DenseMatrix m13(gg, gg_l ,progression_length, transpose);

  mtk::DenseMatrix m14;

  m14 = m13;

  m13.Transpose();

  m14 = m13;

  mtk::Tools::EndUnitTestNo(9);
  mtk::Tools::Assert(m13 == m14);
}

int main () {

  std::cout << "Testing mtk::DenseMatrix class." << std::endl;

  TestDefaultConstructor();
  TestConstructorWithNumRowsNumCols();
  TestConstructAsIdentity();
  TestConstructAsVandermonde();
  TestSetValueGetValue();
  TestConstructAsVandermondeTranspose();
  TestKron();
  TestConstructWithNumRowsNumColsAssignmentOperator();
  TestConstructAsVandermondeTransposeAssignmentOperator();
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
