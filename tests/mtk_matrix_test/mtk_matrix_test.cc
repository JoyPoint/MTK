/*!
\file mtk_matrix_test.cc

\brief Unit test file for the mtk::Matrix class.

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

#include <cmath>
#include <ctime>

#include <iostream>

#include "mtk.h"

void TestDefaultConstructorGetters() {

	mtk::Tools::BeginUnitTestNo(1);

	mtk::Matrix mm;	// A matrix object containing no data at all.

	bool assertion{(mm.encoded_operator() == mtk::EncodedOperator::NOOP) &&
			(mm.storage() == mtk::MatrixStorage::DENSE) &&
			(mm.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
			(mm.num_rows() == 0) &&
			(mm.num_cols() == 0) &&
			(mm.num_values() == 0) &&
			(mm.leading_dimension() == 0) &&
			(mm.num_zero() == 0) &&
			(mm.num_non_zero() == 0) &&
			(mm.num_null() == 0) &&
			(mm.num_non_null() == 0) &&
			(mm.num_low_diags() == 0) &&
			(mm.num_upp_diags() == 0) &&
			(mm.bandwidth() == 0) &&
			(mm.abs_density() == mtk::kZero) &&
			(mm.rel_density() == mtk::kZero) &&
			(mm.abs_sparsity() == mtk::kZero) &&
			(mm.rel_sparsity() == mtk::kZero)};

	mtk::Tools::Assert(assertion);

	mtk::Tools::EndUnitTestNo(1);
}

void TestSettersGettersCopyConstructor() {

	mtk::Tools::BeginUnitTestNo(2);

	const int num_rows{14};
	const int num_cols{17};
	const int num_values{num_rows*num_cols};
	const int num_low_diags{1};
	const int num_upp_diags{1};
	const int bandwidth{3};
	const int num_non_null{};
	const int num_zero{};

	mtk::Matrix m1;	// A matrix object containing no data at all.

	m1.set_encoded_operator(mtk::EncodedOperator::GRADIENT);
	m1.set_storage(mtk::MatrixStorage::DENSE);
	m1.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
	m1.set_num_rows(num_rows);
	m1.set_num_cols(num_cols);
	m1.set_num_low_diags(1);
	m1.set_num_upp_diags(1);

	bool assertion{(m1.encoded_operator() == mtk::EncodedOperator::GRADIENT) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
	  (m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
	  (m1.leading_dimension() == num_cols) &&
	  (m1.num_low_diags() == num_low_diags) &&
	  (m1.num_upp_diags() == num_upp_diags) &&
	  (m1.bandwidth() == bandwidth) &&
	  (m1.num_null() == num_values) &&
	  (m1.num_non_null() == num_non_null) &&
		(m1.abs_density() == mtk::kZero) &&
		(m1.abs_sparsity() == mtk::kOne) &&
	  (m1.num_zero() == num_zero) &&
	  (m1.num_non_zero() == num_values) &&
		(m1.rel_density() == mtk::kZero) &&
		(m1.rel_sparsity() == mtk::kOne)};

	mtk::Matrix m2(m1);	// A matrix object that should be a copy of m1.

	assertion = assertion &&
		(m1.encoded_operator() == m2.encoded_operator()) &&
		(m1.storage() == m2.storage()) &&
		(m1.ordering() == m2.ordering()) &&
		(m1.num_rows() == m2.num_rows()) &&
		(m1.num_cols() == m2.num_cols()) &&
		(m1.num_values() == m2.num_values()) &&
		(m1.leading_dimension() == m2.leading_dimension()) &&
		(m1.num_null() == m2.num_null()) &&
		(m1.num_non_null() == m2.num_non_null()) &&
		(m1.abs_density() == m2.abs_density()) &&
		(m1.abs_sparsity() == m2.abs_sparsity()) &&
		(m1.num_zero() == m2.num_zero()) &&
		(m1.num_non_zero() == m2.num_non_zero()) &&
		(m1.rel_density() == m2.abs_density()) &&
		(m1.rel_sparsity() == m2.rel_sparsity());

	mtk::Tools::Assert(assertion);

	mtk::Tools::EndUnitTestNo(2);
}

void TestIncreaseDecreaseNumNull() {

	mtk::Tools::BeginUnitTestNo(3);

	const int num_rows{19};
	const int num_cols{27};
	const int num_values{num_rows*num_cols};
	const int num_null{num_values - 1};
	const int num_non_null{1};

	mtk::Matrix m1;	// A matrix object containing no data at all.

	mtk::Real abs_density{0.0019493177};
	mtk::Real abs_sparsity{0.99805068};

	m1.set_encoded_operator(mtk::EncodedOperator::DIVERGENCE);
	m1.set_storage(mtk::MatrixStorage::DENSE);
	m1.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
	m1.set_num_rows(num_rows);
	m1.set_num_cols(num_cols);

	bool assertion{(m1.encoded_operator() == mtk::EncodedOperator::DIVERGENCE) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
	  (m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
	  (m1.leading_dimension() == num_cols) &&
	  (m1.num_null() == num_values) &&
	  (m1.num_non_null() == 0) &&
		(m1.abs_density() == mtk::kZero) &&
		(m1.abs_sparsity() == mtk::kOne) &&
	  (m1.num_zero() == 0) &&
	  (m1.num_non_zero() == num_values) &&
		(m1.rel_density() == mtk::kZero) &&
		(m1.rel_sparsity() == mtk::kOne)};

	m1.DecreaseNumNull();

	assertion = assertion &&
		(m1.encoded_operator() == mtk::EncodedOperator::DIVERGENCE) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
		(m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
		(m1.leading_dimension() == num_cols) &&
		(m1.num_null() == num_null) &&
		(m1.num_non_null() == num_non_null) &&
		(std::abs(m1.abs_density() - abs_density) < mtk::kDefaultTolerance) &&
	  (std::abs(m1.abs_sparsity() - abs_sparsity) < mtk::kDefaultTolerance) &&
		(m1.num_zero() == 0) &&
		(m1.num_non_zero() == num_values) &&
		(m1.rel_density() == mtk::kZero) &&
		(m1.rel_sparsity() == mtk::kOne);

	m1.IncreaseNumNull();

  assertion = assertion &&
  	(m1.encoded_operator() == mtk::EncodedOperator::DIVERGENCE) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
	  (m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
	  (m1.leading_dimension() == num_cols) &&
	  (m1.num_null() == num_values) &&
	  (m1.num_non_null() == 0) &&
		(m1.abs_density() == mtk::kZero) &&
		(m1.abs_sparsity() == mtk::kOne) &&
	  (m1.num_zero() == 0) &&
	  (m1.num_non_zero() == num_values) &&
		(m1.rel_density() == mtk::kZero) &&
		(m1.rel_sparsity() == mtk::kOne);

	mtk::Tools::Assert(assertion);

	mtk::Tools::EndUnitTestNo(3);
}

void TestIncreaseDecreaseNumZero() {

	mtk::Tools::BeginUnitTestNo(4);

	const int num_rows{19};
	const int num_cols{27};
	const int num_values{num_rows*num_cols};

	mtk::Matrix m1;	// A matrix object containing no data at all.

	mtk::Real abs_density{0.0019493177};
	mtk::Real abs_sparsity{0.99805068};
	mtk::Real rel_density{0.0019493177};
	mtk::Real rel_sparsity{0.99805068};

	m1.set_encoded_operator(mtk::EncodedOperator::DIVERGENCE);
	m1.set_storage(mtk::MatrixStorage::DENSE);
	m1.set_ordering(mtk::MatrixOrdering::ROW_MAJOR);
	m1.set_num_rows(19);
	m1.set_num_cols(27);

	bool assertion{(m1.encoded_operator() == mtk::EncodedOperator::DIVERGENCE) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
	  (m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
	  (m1.leading_dimension() == num_cols) &&
	  (m1.num_null() == num_values) &&
	  (m1.num_non_null() == 0) &&
		(m1.abs_density() == mtk::kZero) &&
		(m1.abs_sparsity() == mtk::kOne) &&
	  (m1.num_zero() == 0) &&
	  (m1.num_non_zero() == num_values) &&
		(m1.rel_density() == mtk::kZero) &&
		(m1.rel_sparsity() == mtk::kOne)};

	m1.IncreaseNumZero();

	assertion = assertion &&
		(m1.encoded_operator() == mtk::EncodedOperator::DIVERGENCE) &&
		(m1.storage() == mtk::MatrixStorage::DENSE) &&
		(m1.ordering() == mtk::MatrixOrdering::ROW_MAJOR) &&
	  (m1.num_rows() == num_rows) &&
		(m1.num_cols() == num_cols) &&
		(m1.num_values() == num_values) &&
	  (m1.leading_dimension() == num_cols) &&
	  (m1.num_null() == num_values - 1) &&
	  (m1.num_non_null() == 1) &&
		(std::abs(m1.abs_density() - abs_density) < mtk::kDefaultTolerance) &&
	  (std::abs(m1.abs_sparsity() - abs_sparsity) < mtk::kDefaultTolerance) &&
	  (m1.num_zero() == 1) &&
	  (m1.num_non_zero() == num_values - 1) &&
		(std::abs(m1.rel_density() - rel_density) < mtk::kDefaultTolerance) &&
	  (std::abs(m1.rel_sparsity() - rel_sparsity) < mtk::kDefaultTolerance);

	m1.DecreaseNumZero();

	mtk::Tools::Assert(assertion);

	mtk::Tools::EndUnitTestNo(4);
}

int main (int argc, char *argv[]) {

  std::cout << "Testing mtk::Matrix class." << std::endl;

  TestDefaultConstructorGetters();
  TestSettersGettersCopyConstructor();
  TestIncreaseDecreaseNumNull();
  TestIncreaseDecreaseNumZero();
}
