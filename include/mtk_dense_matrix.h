/*!
\file mtk_dense_matrix.h

\brief Defines a common dense matrix, using a 1D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\note We prefer composition to inheritance [Reedy, 2011]. The main reason for
this preference is that inheritance produces a more tightly coupled design. When
a class inherits from another type be it public, protected, or private
inheritance the subclass gains access to all public and protected members of the
base class, whereas with composition, the class is only coupled to the public
members of the other class. Furthermore, if you only hold a pointer to the other
object, then your interface can use a forward declaration of the class rather
than #include its full definition. This results in greater compile-time
insulation and improves the time it takes to compile your code.
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

#ifndef MTK_INCLUDE_DENSE_MATRIX_H_
#define MTK_INCLUDE_DENSE_MATRIX_H_

#include <iostream>

#include "mtk_roots.h"
#include "mtk_enums.h"
#include "mtk_matrix.h"

namespace mtk {

/*!
\class DenseMatrix

\ingroup c04-data_structures

\brief Defines a common dense matrix, using a 1D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.
*/
class DenseMatrix {
 public:
  /// \brief Prints the matrix as a block of numbers (standard way).
  friend std::ostream& operator <<(std::ostream &stream, DenseMatrix &in);

  /// \brief Compares two matrices.
  bool operator ==(const DenseMatrix &in);

  /// \brief Overloaded assignment operator.
  DenseMatrix& operator =(const DenseMatrix &in);

  /// \brief Default constructor.
  DenseMatrix();

  /*!
  \brief Copy constructor.

  \param [in] in Given matrix.
  */
  DenseMatrix(const DenseMatrix &in);

  /*!
  \brief Construct a dense matrix based on the given dimensions.

  \param[in] num_rows Number of rows of the required matrix.
  \param[in] num_cols Number of rows of the required matrix.

  \exception std::bad_alloc
  */
  DenseMatrix(const int &num_rows, const int &num_cols);

  /*!
  \brief Construct a zero-rows-padded identity matrix.

  Used in the construction of the mimetic operators.

  **Def**. A **padded matrix** is a matrix with its first and last rows
  initialized to only zero values:

  \f[
    \bar{\mathbf{I}} = \left(\begin{array}{ccccc}
      0 & 0 & 0 & \dots & 0 \\
      1 & 0 & 0 & \dots & 0 \\
      0 & 1 & 0 & \dots & 0 \\
      \vdots & \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & 0 & \dots & 1 \\
      0 & 0 & 0 & \dots & 0
    \end{array}\right)
  \f]

  \param[in] rank Rank or number of rows/cols in square matrix.
  \param[in] padded Should it be padded?
  \param[in] transpose Should I return the transpose of the requested matrix?

  \exception std::bad_alloc
  */
  DenseMatrix(const int &rank, const bool &padded, const bool &transpose);

  /*!
  \brief Construct a dense Vandermonde matrix.

  **Def**. In linear algebra, a **Vandermonde matrix** is a matrix with
  terms of a geometric progression in each row. This progression uses the terms
  of a given **generator vector**:

  \f[
    \mathbf{V} = \left(\begin{array}{ccccc}
      1 & \alpha_1 & \alpha_1^2 & \dots & \alpha_1^{n-1}\\
      1 & \alpha_2 & \alpha_2^2 & \dots & \alpha_2^{n-1}\\
      1 & \alpha_3 & \alpha_3^2 & \dots & \alpha_3^{n-1}\\
      \vdots & \vdots & \vdots & \ddots &\vdots \\
      1 & \alpha_m & \alpha_m^2 & \dots & \alpha_m^{n-1}
    \end{array}\right)
  \f]

  This constructor generates a Vandermonde matrix, as defined above.

  **Obs**. It in important to understand that the generator vectors to be used
  are nothing but a very particular instance of a grid. These are little chunks,
  little samples, if you will, of a grid which is rectangular and uniform. So
  the selected samples, on the mtk::Div1D and mtk::Grad1D, basically represent
  the entire space, the entire grid. This is why nor the CRS nor the CBS
  algorithms may work for irregular geometries, such as curvilinear grids.

  \param[in] gen Given generator vector.
  \param[in] gen_length Length generator vector.
  \param[in] pro_length Length the progression.
  \param[in] transpose Should the transpose be created instead?

  \exception std::bad_alloc
  */
  DenseMatrix(const Real *gen,
              const int &gen_length,
              const int &pro_length,
              const bool &transpose);

  /// \brief Destructor.
  ~DenseMatrix();

  /*!
  \brief Provides access to the matrix data.

  \return Pointer to a Matrix.
  */
  Matrix matrix_properties() const;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_rows() const;

  /*!
  \brief Gets the number of columns.

  \return Number of columns of the matrix.
  */
  int num_cols() const;

  /*!
  \brief Provides access to the matrix value array.

  \return Pointer to an array of mtk::Real.
  */
  Real* data() const;

   /*!
  \brief Sets the ordering of the matrix.

  \param[in] oo Ordering.

  \return The required value at the specified coordinates.
  */
  void SetOrdering(mtk::MatrixOrdering oo);

  /*!
  \brief Gets a value on the given coordinates.

  \param[in] row_coord Row coordinate.
  \param[in] col_coord Column coordinate.

  \return The required value at the specified coordinates.
  */
  Real GetValue(const int &row_coord, const int &col_coord) const;

  /*!
  \brief Sets a value on the given coordinates.

  \param[in] row_coord  Row coordinate.
  \param[in] col_coord  Column coordinate.
  \param[in] val        Row Actual value to be inserted.
  */
  void SetValue(const int &row_coord,
                const int &col_coord,
                const Real &val);

  /// \brief Transpose this matrix.
  void Transpose();

  /// \brief Make the matrix row-wise ordered.
  void OrderRowMajor();

  /// \brief Make the matrix column-wise ordered.
  void OrderColMajor();

  /*!
  \brief Construct a dense matrix based on the Kronecker product of arguments.

  \param[in] aa First matrix.
  \param[in] bb Second matrix.

  \exception std::bad_alloc

  \todo Implement Kronecker product using the BLAS.
  */
  static DenseMatrix Kron(const DenseMatrix &aa, const DenseMatrix &bb);

 private:
  Matrix matrix_properties_;  ///< Data related to the matrix nature.

  Real *data_; ///< Array holding the data in contiguouos position in memory.
};
}
#endif  // End of: MTK_INCLUDE_MTK_DENSE_MATRIX_H_
