/*!
\file mtk_matrix.h

\brief Definition of the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.

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

#ifndef MTK_INCLUDE_MATRIX_H_
#define	MTK_INCLUDE_MATRIX_H_

#include <iostream>

#include "mtk_roots.h"
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
  /// \brief Default constructor.
  Matrix();

  /*!
  \brief Copy constructor.

  \param[in] in Given matrix.
  */
  Matrix(const Matrix &in);

  /// \brief Destructor.
  ~Matrix();

  /// \brief Gets the type of storage of this matrix.
  MatrixStorage storage() const;

  /// \brief Gets the ordering of this matrix.
  MatrixOrdering ordering() const;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_rows() const;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_cols() const;

  /*!
  \brief Gets the number of values.

  \return Number of values of the matrix.
  */
  int num_values() const;

  /*!
  \brief Gets the matrix' leading dimension.

  Leading dimension of the data array is the number of elements between
  successive rows (for row major storage) in memory. Most of the cases, the
  leading dimension is the same as the number of columns.

  \return Leading dimension of the matrix.
  */
  int ld() const;

  /*!
  \brief Gets the number of zeros.

  \return Number of zeros of the matrix.
  */
  int num_zero() const;

  /*!
  \brief Gets the number of non-zero values.

  \return Number of non-zero values of the matrix.
  */
  int num_non_zero() const;

  /*!
  \brief Gets the number of null values.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Number of null values of the matrix.
  */
  int num_null() const;

  /*!
  \brief Gets the number of non-null values.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Number of non-null values of the matrix.
  */
  int num_non_null() const;
  
  /*!
  \brief Gets the number of lower diagonals.

  \return Number of lower diagonals.
  */
  int kl() const;

  /*!
  \brief Gets the number of upper diagonals.

  \return Number of upper diagonals.
  */
  int ku() const;

  /*!
  \brief Gets the bandwidth.

  \return Bandwidth of the matrix.
  */
  int bandwidth() const;

  /*!
  \brief Gets the absolute density.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Absolute density of the matrix.
  */
  Real abs_density() const;

  /*!
  \brief Gets the relative density.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Relative density of the matrix.
  */
  Real rel_density() const;

  /*!
  \brief Gets the Absolute sparsity.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Absolute sparsity of the matrix.
  */
  Real abs_sparsity() const;

  /*!
  \brief Gets the Relative sparsity.

  \sa http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf

  \return Relative sparsity of the matrix.
  */
  Real rel_sparsity() const;

  /*!
  \brief Sets the storage type of the matrix.

  \sa MatrixStorage

  \param[in] tt Type of the matrix storage.
  */
  void set_storage(const MatrixStorage &tt);

  /*!
  \brief Sets the ordering of the matrix.

  \sa MatrixOrdering

  \param[in] oo Ordering of the matrix.
  */
  void set_ordering(const MatrixOrdering &oo);

  /*!
  \brief Sets the number of rows of the matrix.

  \param[in] num_rows Number of rows.
  */
  void set_num_rows(int num_rows);

  /*!
  \brief Sets the number of columns of the matrix.

  \param[in] num_cols Number of columns.
  */
  void set_num_cols(int num_cols);

  /*!
  \brief Sets the number of zero values of the matrix that matter.

  \param[in] in Number of zero values.
  */
  void set_num_zero(int in);

  /*!
  \brief Sets the number of zero values of the matrix that DO NOT matter.

  \param[in] in Number of zero values.
  */
  void set_num_null(int in);

  /// \brief Increases the number of values that equal zero but with meaning.
  void IncreaseNumZero();

  /// \brief Increases the number of values that equal zero but with no meaning.
  void IncreaseNumNull();

 private:
  MatrixStorage storage_; ///< What type of matrix is this?

  MatrixOrdering ordering_; ///< What kind of ordering is it following?

  int num_rows_;        ///< Number of rows.
  int num_cols_;        ///< Number of columns.
  int num_values_;      ///< Number of total values in matrix.
  int ld_;              ///< Elements between successive rows when row-major.

  int num_zero_;        ///< Number of zeros.
  int num_non_zero_;    ///< Number of non-zero values.
  int num_null_;        ///< Number of null (insignificant) values.
  int num_non_null_;    ///< Number of null (significant) values.

  int kl_;              ///< Number of lower diagonals on a banded matrix.
  int ku_;              ///< Number of upper diagonals on a banded matrix.
  int bandwidth_;       ///< Bandwidth of the matrix.

  Real abs_density_;    ///< Absolute density of matrix.
  Real rel_density_;    ///< Relative density of matrix.
  Real abs_sparsity_;   ///< Absolute sparsity of matrix.
  Real rel_sparsity_;   ///< Relative sparsity of matrix.
};
}
#endif  // End of: MTK_INCLUDE_MATRIX_H_
