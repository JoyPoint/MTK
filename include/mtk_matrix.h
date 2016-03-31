/*!
\file mtk_matrix.h

\brief Definition of the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.

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

#ifndef MTK_INCLUDE_MATRIX_H_
#define	MTK_INCLUDE_MATRIX_H_

#include <iostream>

#include "mtk_foundations.h"
#include "mtk_enums.h"

namespace mtk {

/*!
\class Matrix

\ingroup c04-data_structures

\brief Definition of the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.
*/
class Matrix {
 public:
  /*!
	\brief Default constructor.
	*/
  Matrix();

  /*!
  \brief Copy constructor.

  \param[in] in Given matrix.
  */
  Matrix(const Matrix &in);

  /*!
  \brief Destructor.
  */
  ~Matrix() noexcept;

  /*!
  \brief Gets the type of mimetic operator encoded by this matrix.

  \return Type of mimetic operator encoded by this matrix.

  \sa mtk::EncodedOperator.
  */
  EncodedOperator encoded_operator() const noexcept;

  /*!
  \brief Gets the type of storage of this matrix.

  \return Type of storage of this matrix.

  \sa mtk::MatrixStorage.
  */
  MatrixStorage storage() const noexcept;

  /*!
  \brief Gets the type of ordering of this matrix.

  \return Type of ordering of this matrix.

  \sa mtk::MatrixOrdering.
  */
  MatrixOrdering ordering() const noexcept;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_rows() const noexcept;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_cols() const noexcept;

  /*!
  \brief Gets the number of values.

  \return Number of values of the matrix.
  */
  int num_values() const noexcept;

  /*!
  \brief Gets the matrix' leading dimension.

  Leading dimension of the data array is the number of elements between
  successive rows (for row major storage) in memory. Most of the cases, the
  leading dimension is the same as the number of columns.

  \return Leading dimension of the matrix.
  */
  int leading_dimension() const noexcept;

  /*!
  \brief Gets the number of zeros.

  \return Number of zeros of the matrix.
  */
  int num_zero() const noexcept;

  /*!
  \brief Gets the number of non-zero values.

  \return Number of non-zero values of the matrix.
  */
  int num_non_zero() const noexcept;

  /*!
  \brief Gets the number of null values.

  \return Number of null values of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  int num_null() const noexcept;

  /*!
  \brief Gets the number of non-null values.

  \return Number of non-null values of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  int num_non_null() const noexcept;

  /*!
  \brief Gets the number of lower diagonals.

  \return Number of lower diagonals.
  */
  int num_low_diags() const noexcept;

  /*!
  \brief Gets the number of upper diagonals.

  \return Number of upper diagonals.
  */
  int num_upp_diags() const noexcept;

  /*!
  \brief Gets the bandwidth.

  \return Bandwidth of the matrix.
  */
  int bandwidth() const noexcept;

  /*!
  \brief Gets the absolute density.

  \return Absolute density of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  Real abs_density() const noexcept;

  /*!
  \brief Gets the relative density.

  \return Relative density of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  Real rel_density() const noexcept;

  /*!
  \brief Gets the Absolute sparsity.

  \return Absolute sparsity of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  Real abs_sparsity() const noexcept;

  /*!
  \brief Gets the Relative sparsity.

  \return Relative sparsity of the matrix.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
  */
  Real rel_sparsity() const noexcept;

  /*!
  \brief Sets the type of encoded operator of the matrix.

  \param[in] in Type of encoded operator.

  \sa mtk::EncodedOperator
  */
  void set_encoded_operator(const EncodedOperator &in) noexcept;

  /*!
  \brief Sets the storage type of the matrix.

  \param[in] tt Type of the matrix storage.

  \sa mtk::MatrixStorage
  */
  void set_storage(const MatrixStorage &tt) noexcept;

  /*!
  \brief Sets the ordering of the matrix.

  \param[in] oo Ordering of the matrix.

  \sa mtk::MatrixOrdering
  */
  void set_ordering(const MatrixOrdering &oo) noexcept;

  /*!
  \brief Sets the number of rows of the matrix.

  \param[in] num_rows Number of rows.
  */
  void set_num_rows(const int &num_rows) noexcept;

  /*!
  \brief Sets the number of columns of the matrix.

  \param[in] num_cols Number of columns.
  */
  void set_num_cols(const int &num_cols) noexcept;

  /*!
  \brief Sets the number of lower diagonals of the matrix.

  \param[in] num_low_diags Number of lower diagonals.
  */
  void set_num_low_diags(const int &num_low_diags) noexcept;

  /*!
  \brief Sets the number of upper diagonals of the matrix.

  \param[in] num_upp_diags Number of upper diagonals.
  */
  void set_num_upp_diags(const int &num_upp_diags) noexcept;

  /*!
  \brief Sets the number of null elements of the matrix.

  \param[in] num_null Number of null elements.
  */
  void set_num_null(const int &num_null) noexcept;

  /*!
  \brief Sets the number of zeros of the matrix.

  \param[in] num_zero Number of zeros on the matrix.
  */
  void set_num_zero(const int &num_zero) noexcept;

  /*!
  \brief Decreases the number of values that equal zero but with no meaning.
  */
  void IncreaseNumNull() noexcept;

  /*!
	\brief Decreases the number of values that equal zero but with no meaning.
  */
  void DecreaseNumNull() noexcept;

  /*!
  \brief Increases the number of values that equal zero but with meaning.
  */
  void IncreaseNumZero() noexcept;

  /*!
  \brief Decreases the number of values that equal zero but with meaning.
  */
  void DecreaseNumZero() noexcept;

 private:
  /*!
	\brief Computes the leading dimension of the matrix.

	\param[in] num_rows Number of rows.
	\param[in] num_cols Number of columns.

	\return Number of values of the matrix.
	*/
  int ComputeNumValues(const int &num_rows, const int &num_cols) const noexcept;

  /*!
  \brief Computes the leading dimension of the matrix.

	\param[in] num_rows Number of rows.
  \param[in] num_cols Number of columns.

  \return Leading dimension of the matrix.
  */
  int ComputeLeadingDimension(const int &num_rows, const int &num_cols)
  	const noexcept;

  /*!
	\brief Computes the bandwidth of the matrix.

	\param[in] num_low_diags Number of lower diagonals.
	\param[in] num_upp_diags Number of upper diagonals.

	\return Bandwidth of the matrix.
	*/
  int ComputeBandwidth(const int &num_low_diags, const int &num_upp_diags)
  	const noexcept;

  /*!
	\brief Computes the number of non-null values of the matrix.

	\param[in] num_values Number of values of the matrix.
	\param[in] num_null Number of null values of the matrix.

	\return Number of non-null values of the matrix.
	*/
  int ComputeNumNonNull(const int &num_values, const int &num_null)
  	const noexcept;

  /*!
	\brief Computes the absolute density of the matrix.

	Defined as
	\f[
		\frac{\textrm{num_non_null}}{\textrm{num_values}}.
	\f]

	\param[in] num_non_null Number of non-null values of the matrix.
	\param[in] num_values Number of total values of the matrix.

	\return Absolute density of the matrix.
	*/
  mtk::Real ComputeAbsDensity(const int &num_non_null, const int &num_values)
  	const noexcept;

  /*!
	\brief Computes the absolute sparsity of the matrix.

	Defined as
	\f[
		1 - \frac{\textrm{num_non_null}}{\textrm{num_values}}.
	\f]

	\param[in] absolute_density Absolute density of the matrix.

	\return Absolute sparsity of the matrix.
	*/
  mtk::Real ComputeAbsSparsity(const mtk::Real &absolute_density)
  	const noexcept;

  /*!
	\brief Computes the number of non-zero values of the matrix.

	\param[in] num_values Number of values of the matrix.
	\param[in] num_zero Number of zero values of the matrix.

	\return Number of non-zero values of the matrix.
	*/
  int ComputeNumNonZero(const int &num_values, const int &num_zero)
  	const noexcept;

  /*!
	\brief Computes the relative density of the matrix.

	Defined as
	\f[
		\frac{\textrm{num_non_zero}}{\textrm{num_values}}.
	\f]

	\param[in] num_non_zero Number of non-zero values of the matrix.
	\param[in] num_values Number of total values of the matrix.

	\return Relative density of the matrix.
	*/
  mtk::Real ComputeRelDensity(const int &num_non_zero, const int &num_values)
  	const noexcept;

  /*!
	\brief Computes the relative sparsity of the matrix.

	Defined as
	\f[
		1 - \frac{\textrm{num_non_zero}}{\textrm{num_values}}.
	\f]

	\param[in] relative_density Relative density of the matrix.

	\return Relative sparsity of the matrix.
	*/
  mtk::Real ComputeRelSparsity(const mtk::Real &relative_density)
  	const noexcept;

  EncodedOperator encoded_operator_;	/// Type of mimetic operator encoded.

  MatrixStorage storage_; ///< What type of matrix is this?

  MatrixOrdering ordering_; ///< What kind of ordering is it following?

  int num_rows_;        	///< Number of rows.
  int num_cols_;        	///< Number of columns.
  int num_values_;      	///< Number of total values in matrix.
  int leading_dimension_;	///< Elements between successive rows when row-major.

  int num_low_diags_;	///< Number of lower diagonals on a banded matrix.
  int num_upp_diags_; ///< Number of upper diagonals on a banded matrix.
  int bandwidth_;     ///< Bandwidth of the matrix.

  int num_null_;      ///< Number of null (no meaning) values.
  int num_non_null_;	///< Number of null values.

  Real abs_density_;  ///< Absolute density of matrix.
  Real abs_sparsity_; ///< Absolute sparsity of matrix.

  int num_zero_;      ///< Number of zeros (with meaning).
  int num_non_zero_;  ///< Number of non-zero values.

  Real rel_density_;  ///< Relative density of matrix.
  Real rel_sparsity_;	///< Relative sparsity of matrix.
};
}
#endif  // End of: MTK_INCLUDE_MATRIX_H_
