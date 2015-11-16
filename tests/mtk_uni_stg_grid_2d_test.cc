/*!
\file mtk_uni_stg_grid_2d_test.cc

\brief Test file for the mtk::UniStgGrid2D class.

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

#include <cmath>
#include <ctime>

#include <iostream>

#include "mtk.h"

void TestDefaultConstructor() {

  mtk::Tools::BeginUnitTestNo(1);

  mtk::UniStgGrid2D gg;

  mtk::Tools::EndUnitTestNo(1);
  mtk::Tools::Assert(gg.delta_x() == mtk::kZero && gg.delta_y() == mtk::kZero);
}

void
TestConstructWithWestEastNumCellsXSouthNorthBndysNumCellsYOStreamOperator() {

  mtk::Tools::BeginUnitTestNo(2);

  mtk::Real aa = 0.0;
  mtk::Real bb = 1.0;
  mtk::Real cc = 0.0;
  mtk::Real dd = 1.0;

  int nn = 5;
  int mm = 7;

  mtk::UniStgGrid2D gg(aa, bb, nn, cc, dd, mm);

  std::cout << gg << std::endl;

  mtk::Tools::EndUnitTestNo(2);
  mtk::Tools::Assert(gg.delta_x() == 0.2 &&
                     abs(gg.delta_y() - 0.142857) < mtk::kDefaultTolerance);
}

void TestGetters() {

  mtk::Tools::BeginUnitTestNo(3);

  mtk::Real aa = 0.0;
  mtk::Real bb = 1.0;
  mtk::Real cc = 0.0;
  mtk::Real dd = 1.0;

  int nn = 5;
  int mm = 7;

  mtk::UniStgGrid2D gg(aa, bb, nn, cc, dd, mm);

  bool assertion{true};

  assertion = assertion && (gg.west_bndy() == aa);
  assertion = assertion && (gg.east_bndy() == bb);
  assertion = assertion && (gg.num_cells_x() == nn);
  assertion = assertion && (gg.south_bndy() == cc);
  assertion = assertion && (gg.north_bndy() == dd);
  assertion = assertion && (gg.num_cells_y() == mm);

  mtk::Tools::EndUnitTestNo(3);
  mtk::Tools::Assert(assertion);
}

mtk::Real ScalarField(mtk::Real xx, mtk::Real yy) {

  mtk::Real aux{-(1.0/2.0)*xx*xx - (1.0/2.0)*yy*yy};

  return xx*yy*exp(aux);
}

void TestBindScalarFieldWriteToFile() {

  mtk::Tools::BeginUnitTestNo(4);

  mtk::Real aa = 0.0;
  mtk::Real bb = 1.0;
  mtk::Real cc = 0.0;
  mtk::Real dd = 1.0;

  int nn = 5;
  int mm = 5;

  mtk::UniStgGrid2D gg(aa, bb, nn, cc, dd, mm);

  gg.BindScalarField(ScalarField);

  if(!gg.WriteToFile("mtk_uni_stg_grid_2d_test_04.dat", "x", "y", "u(x,y)")) {
    std::cerr << "Error writing to file." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(4);
}

mtk::Real VectorFieldPComponent(mtk::Real xx, mtk::Real yy) {

  return xx + 0.01;
}

mtk::Real VectorFieldQComponent(mtk::Real xx, mtk::Real yy) {

  return yy + 0.01;
}

void TestBindVectorField() {

  mtk::Tools::BeginUnitTestNo(5);

  mtk::Real aa = 0.0;
  mtk::Real bb = 1.0;
  mtk::Real cc = 0.0;
  mtk::Real dd = 1.0;

  int nn = 5;
  int mm = 5;

  mtk::UniStgGrid2D gg(aa, bb, nn, cc, dd, mm, mtk::VECTOR);

  gg.BindVectorField(VectorFieldPComponent, VectorFieldQComponent);

  std::cout << gg << std::endl;

  if(!gg.WriteToFile("mtk_uni_stg_grid_2d_test_05.dat", "x", "y", "v(x,y)")) {
    std::cerr << "Error writing to file." << std::endl;
  }

  mtk::Tools::EndUnitTestNo(5);
}

int main () {

  std::cout << "Testing mtk::UniStgGrid2D class." << std::endl;

  TestDefaultConstructor();
  TestConstructWithWestEastNumCellsXSouthNorthBndysNumCellsYOStreamOperator();
  TestGetters();
  TestBindScalarFieldWriteToFile();
  TestBindVectorField();
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
