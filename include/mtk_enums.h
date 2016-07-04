/*!
\file mtk_enums.h

\brief Definitions of the enumeration types in the MTK.

Enumeration types are used throughout the MTK to differentiate instances of
derived classes, as well as for mnemonic purposes. In this file, definitions for
the enumeration types are listed alphabetically.

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

#ifndef MTK_INCLUDE_ENUMS_H_
#define MTK_INCLUDE_ENUMS_H_

namespace mtk {

/*!
\enum MatrixStorage

\ingroup c02-enums

\brief Considered matrix storage schemes to implement sparse matrices.

The considered sparse storage schemes are selected so that these are compatible
with some of the most used mathematical APIs, as follows: DENSE and BANDED for
<a href="http://www.netlib.org/blas/">BLAS</a>,
<a href="http://www.netlib.org/lapack/">LAPACK</a>, and
<a href="http://www.netlib.org/scalapack/">ScaLAPACK</a>. Finally, CRS for
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU</a>.
*/
enum class MatrixStorage {
  DENSE,    ///< Dense matrices, implemented as a 1D array: DenseMatrix.
  BANDED,   ///< Banded matrices ala LAPACK and ScaLAPACK: Must be implemented.
  CRS       ///< Compressed-Rows Storage: Must be implemented.
};

/*!
\enum MatrixOrdering

\ingroup c02-enums

\brief Considered matrix ordering (for Fortran purposes).

Row-major ordering is used for most application in C/C++. For Fortran purposes,
the matrices must be listed in a column-major ordering.

\sa https://en.wikipedia.org/wiki/Row-major_order
*/
enum class MatrixOrdering {
  ROW_MAJOR,  ///< Row-major ordering (C/C++).
  COL_MAJOR   ///< Column-major ordering (Fortran).
};

/*!
\enum FieldNature

\ingroup c02-enums

\brief Nature of the field discretized in a given grid.

Fields can be **scalar** or **vector** in nature.

\sa https://en.wikipedia.org/wiki/Scalar_field

\sa https://en.wikipedia.org/wiki/Vector_field
*/
enum class FieldNature {
  SCALAR,  ///< Scalar-valued field.
  VECTOR   ///< Vector-valued field.
};

/*!
\enum DirInterp

\ingroup c02-enums

\brief Interpolation operator.

Used to tag different directions of interpolation supported.
*/
enum class DirInterp {
  SCALAR_TO_VECTOR, ///< Interpolations places scalar on vectors' location.
  VECTOR_TO_SCALAR  ///< Interpolations places vectors on scalars' location.
};

/*!
\enum EncodedOperator

\ingroup c02-enums

\brief Operators matrices can encode.

Used to tag different different operators a matrix can encode.
*/
enum class EncodedOperator {
  NOOP,             ///< No operator.
  GRADIENT,         ///< Gradient operator.
  DIVERGENCE,       ///< Divergence operator.
  INTERPOLATION,    ///< Interpolation operator.
  CURL,             ///< Curl operator.
  LAPLACIAN         ///< Laplacian operator.
};
}
#endif  // End of: MTK_INCLUDE_ENUMS_H_
