/*!
\file mtk_crs_matrix_test.cc

\brief Test file for the mtk::CRSMatrix class.

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

#include <vector>

#include "mtk_tools.h"
#include "mtk_crs_matrix.h"

void TestMatrixVectorMultiplication1() {

	mtk::Tools::BeginUnitTestNo(1);

	mtk::CRSMatrix m1(3);

	m1.SetValue(1.0, 0, 2);
	m1.SetValue(2.0, 1, 0);
	m1.SetValue(3.0, 1, 2);
	m1.SetValue(4.0, 2, 1);
	m1.SetValue(7.0, 2, 2);

	std::vector<mtk::Real> x1{};

	x1.push_back(2.0);
	x1.push_back(0.0);
	x1.push_back(3.0);

	std::vector<mtk::Real> result1{};

	result1.push_back(3.0);
	result1.push_back(13.0);
	result1.push_back(21.0);

	std::vector<mtk::Real> rr{};

	rr = m1.Multiply(x1);

	mtk::Tools::EndUnitTestNo(1);
	mtk::Tools::Assert(rr == result1);
}

void TestMatrixVectorMultiplication2() {

	mtk::Tools::BeginUnitTestNo(2);

	mtk::CRSMatrix m2(3);
	std::vector<mtk::Real> x2(3, 1);
	std::vector<mtk::Real> result2(3, 0);

	mtk::Tools::EndUnitTestNo(2);
	mtk::Tools::Assert(m2.Multiply(x2) == result2);
}

void TestMatrixVectorMultiplication3() {

	mtk::Tools::BeginUnitTestNo(3);

	mtk::CRSMatrix m3(3);
	m3.SetValue(1.0, 0, 0)
		.SetValue(2.0, 0, 1)
		.SetValue(3.0, 0, 2)
		.SetValue(4.0, 1, 0)
		.SetValue(5.0, 1, 1)
		.SetValue(6.0, 1, 2)
		.SetValue(7.0, 2, 0)
		.SetValue(8.0, 2, 1)
		.SetValue(9.0, 2, 2);

	std::vector<mtk::Real> x3(3, 0.0);
	std::vector<mtk::Real> result3(3, 0.0);

	mtk::Tools::EndUnitTestNo(3);
	mtk::Tools::Assert(m3.Multiply(x3) == result3);
}

void TestMatrixMatrixAdition1() {

	mtk::Tools::BeginUnitTestNo(4);

	mtk::CRSMatrix a4(3);
	mtk::CRSMatrix b4(3);
	b4.SetValue(1.0, 0, 0)
		.SetValue(2.0, 0, 1)
		.SetValue(3.0, 0, 2)
		.SetValue(4.0, 1, 0)
		.SetValue(5.0, 1, 1)
		.SetValue(6.0, 1, 2)
		.SetValue(7.0, 2, 0)
		.SetValue(8.0, 2, 1)
		.SetValue(9.0, 2, 2);

	mtk::Tools::EndUnitTestNo(4);
	mtk::Tools::Assert(a4.Add(b4) == b4);
}

void TestMatrixMatrixAdition2() {

	mtk::Tools::BeginUnitTestNo(5);

	mtk::CRSMatrix a5(3);
	a5.SetValue(-1.0, 0, 0)
		.SetValue(-2.0, 0, 1)
		.SetValue(-3.0, 0, 2)
		.SetValue(-4.0, 1, 0)
		.SetValue(-5.0, 1, 1)
		.SetValue(-6.0, 1, 2)
		.SetValue(-7.0, 2, 0)
		.SetValue(-8.0, 2, 1)
		.SetValue(-9.0, 2, 2);

	mtk::CRSMatrix b5(3);
	b5.SetValue(1.0, 0, 0)
		.SetValue(2.0, 0, 1)
		.SetValue(3.0, 0, 2)
		.SetValue(4.0, 1, 0)
		.SetValue(5.0, 1, 1)
		.SetValue(6.0, 1, 2)
		.SetValue(7.0, 2, 0)
		.SetValue(8.0, 2, 1)
		.SetValue(9.0, 2, 2);

	mtk::Tools::EndUnitTestNo(5);
	mtk::Tools::Assert(a5.Add(b5) == mtk::CRSMatrix(3));
}

void TestMatrixMatrixAdition3() {

	mtk::Tools::BeginUnitTestNo(6);

	mtk::CRSMatrix a6(3);
	a6.SetValue(1.0, 0, 0)
		.SetValue(1.0, 0, 2)
		.SetValue(1.0, 1, 1)
		.SetValue(1.0, 2, 0)
		.SetValue(1.0, 2, 2);

	mtk::CRSMatrix b6(3);
	b6.SetValue(1.0, 0, 1)
		.SetValue(1.0, 1, 0)
		.SetValue(1.0, 1, 2)
		.SetValue(1.0, 2, 1);

	mtk::CRSMatrix result6(3);
	result6.SetValue(1.0, 0, 0)
				 .SetValue(1.0, 0, 1)
				 .SetValue(1.0, 0, 2)
				 .SetValue(1.0, 1, 0)
				 .SetValue(1.0, 1, 1)
				 .SetValue(1.0, 1, 2)
				 .SetValue(1.0, 2, 0)
				 .SetValue(1.0, 2, 1)
				 .SetValue(1.0, 2, 2);

	mtk::Tools::EndUnitTestNo(6);
	mtk::Tools::Assert(a6.Add(b6) == result6);
}

void TestMatrixMatrixMultiplication1() {

	mtk::Tools::BeginUnitTestNo(7);

	mtk::CRSMatrix a7(3);
	a7.SetValue(1.0, 0, 0)
		.SetValue(2.0, 0, 1)
		.SetValue(3.0, 0, 2)
		.SetValue(4.0, 1, 0)
		.SetValue(5.0, 1, 1)
		.SetValue(6.0, 1, 2)
		.SetValue(7.0, 2, 0)
		.SetValue(8.0, 2, 1)
		.SetValue(9.0, 2, 2);

	mtk::CRSMatrix b7(3);
	b7.SetValue(1.0, 0, 0)
		.SetValue(1.0, 0, 1)
		.SetValue(1.0, 0, 2)
		.SetValue(1.0, 1, 0)
		.SetValue(1.0, 1, 1)
		.SetValue(1.0, 1, 2)
		.SetValue(1.0, 2, 1);

	mtk::CRSMatrix result7(3);
	result7.SetValue(3.0, 0, 0)
				 .SetValue(6.0, 0, 1)
				 .SetValue(3.0, 0, 2)
				 .SetValue(9.0, 1, 0)
				 .SetValue(15.0, 1, 1)
				 .SetValue(9.0, 1, 2)
				 .SetValue(15.0, 2, 0)
				 .SetValue(24.0, 2, 1)
				 .SetValue(15.0, 2, 2);

	mtk::Tools::EndUnitTestNo(7);
	mtk::Tools::Assert(a7.Multiply(b7) == result7);
}

void TestMatrixMatrixMultiplication2() {

	mtk::Tools::BeginUnitTestNo(8);

	mtk::CRSMatrix a8(3);
	a8.SetValue(1.0, 0, 0)
		.SetValue(1.0, 0, 2)
		.SetValue(1.0, 1, 2)
		.SetValue(2.0, 2, 0)
		.SetValue(2.0, 2, 1);

	mtk::CRSMatrix b8(3);
	mtk::CRSMatrix result8(3);

	mtk::Tools::EndUnitTestNo(8);
	mtk::Tools::Assert(a8.Multiply(b8) == result8);
}

void TestMatrixMatrixMultiplication3() {

	mtk::Tools::BeginUnitTestNo(9);

	mtk::CRSMatrix a9(3);
	a9.SetValue(1.0, 0, 0)
		.SetValue(1.0, 1, 1)
		.SetValue(1.0, 2, 2);

	mtk::CRSMatrix b9(3);
	b9.SetValue(1.0, 0, 0)
		.SetValue(2.0, 0, 1)
		.SetValue(3.0, 0, 2)
		.SetValue(4.0, 1, 0)
		.SetValue(5.0, 1, 1)
		.SetValue(6.0, 1, 2)
		.SetValue(7.0, 2, 0)
		.SetValue(8.0, 2, 1)
		.SetValue(9.0, 2, 2);

	mtk::CRSMatrix result9(3);
	result9.SetValue(1.0, 0, 0)
				 .SetValue(2.0, 0, 1)
				 .SetValue(3.0, 0, 2)
				 .SetValue(4.0, 1, 0)
				 .SetValue(5.0, 1, 1)
				 .SetValue(6.0, 1, 2)
				 .SetValue(7.0, 2, 0)
				 .SetValue(8.0, 2, 1)
				 .SetValue(9.0, 2, 2);

	mtk::Tools::EndUnitTestNo(9);
	mtk::Tools::Assert(a9.Multiply(b9) == result9 && b9.Multiply(a9) == result9);
}

void TestKron() {

	mtk::Tools::BeginUnitTestNo(10);

  mtk::CRSMatrix a10(6);

  a10.SetValue(10.0, 0, 0)
     .SetValue(-2.0, 0, 4)
     .SetValue(3.0, 1, 0)
     .SetValue(9.0, 1, 1)
     .SetValue(3.0, 1, 5)
     .SetValue(7.0, 2, 1)
     .SetValue(8.0, 2, 4)
     .SetValue(7.0, 2, 3)
     .SetValue(3.0, 3, 0)
     .SetValue(8.0, 3, 2)
     .SetValue(7.0, 3, 3)
     .SetValue(5.0, 3, 4)
     .SetValue(8.0, 4, 1)
     .SetValue(9.0, 4, 3)
     .SetValue(9.0, 4, 4)
     .SetValue(13.0, 4, 5)
     .SetValue(4.0, 5, 1)
     .SetValue(2.0, 5, 4)
     .SetValue(-1.0, 5, 5);

  mtk::CRSMatrix a11(3);
  a11.SetValue(1.0, 0, 0)
     .SetValue(1.0, 1, 1)
     .SetValue(1.0, 2, 2);

  mtk::CRSMatrix b11(2);

  b11.SetValue(2.0, 0, 0)
     .SetValue(2.0, 0, 1)
     .SetValue(2.0, 1, 0)
     .SetValue(2.0, 1, 1);

  mtk::CRSMatrix result11(6);

  result11.SetValue(2.0, 0, 0)
    		  .SetValue(2.0, 0, 1)
    		  .SetValue(2.0, 1, 0)
    		  .SetValue(2.0, 1, 1)
    		  .SetValue(2.0, 2, 2)
    		  .SetValue(2.0, 2, 3)
    		  .SetValue(2.0, 3, 2)
    		  .SetValue(2.0, 3, 3)
    		  .SetValue(2.0, 4, 4)
    		  .SetValue(2.0, 4, 5)
    		  .SetValue(2.0, 5, 4)
    		  .SetValue(2.0, 5, 5);

	mtk::Tools::EndUnitTestNo(10);
	mtk::Tools::Assert(mtk::CRSMatrix::Kron(a11, b11) == result11);

}

int main(int argc, char ** argv) {

	std::cout << "Testing mtk::CRSMatrix class." << std::endl;

	TestMatrixVectorMultiplication1();
	TestMatrixVectorMultiplication2();
	TestMatrixVectorMultiplication3();
	TestMatrixMatrixAdition1();
	TestMatrixMatrixAdition2();
	TestMatrixMatrixAdition3();
	TestMatrixMatrixMultiplication1();
	TestMatrixMatrixMultiplication2();
	TestMatrixMatrixMultiplication3();
	TestKron();
}
