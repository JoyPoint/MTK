/*!
\file mtk_blas_facade.cc

\brief Adapter class for the BLAS API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the BLAS.

The **BLAS (Basic Linear Algebra Subprograms)** are routines that provide
standard building blocks for performing basic vector and matrix operations. The
Level 1 BLAS perform scalar, vector and vector-vector operations, the Level 2
BLAS perform matrix-vector operations, and the Level 3 BLAS perform
matrix-matrix operations.

The BLAS can be installed from links given in the See Also section of this page.

\sa http://www.netlib.org/blas/

\sa https://software.intel.com/en-us/non-commercial-software-development

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

#include <iostream>
#include <iomanip>

#include <vector>

#include "mtk_roots.h"
#include "mtk_tools.h"
#include "mtk_blas_adapter.h"

namespace mtk {

extern "C" {

#ifdef MTK_PRECISION_DOUBLE
/*!
\fn dnrm2_

\brief Double-precision computation of the norm 2:

\sa http://www.netlib.org/lapack/explore-html/da/d7f/dnrm2_8f.html

\param[in] n    Length of the input array.
\param[in] x    Given array.
\param[in] incx Increment of the values in the given array.

\return Double-precision computation of the norm 2.
*/
double dnrm2_(int *n, double *x, int *incx);
#else
/*!
\fn snrm2_

\brief Single-precision computation of the norm 2:

\sa http://www.netlib.org/lapack/explore-html/d7/df1/snrm2_8f.html

\param[in] n    Length of the input array.
\param[in] x    Given array.
\param[in] incx Increment of the values in the given array.

\return Single-precision computation of the norm 2.
*/
float snrm2_(int *n, float *x, int *incx);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\fn daxpy_

\brief DAXPY constant times a vector plus a vector.

DAXPY constant times a vector plus a vector. Uses unrolled loops for increments
equal to one.

\sa http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html

\param[in] n    Length of the input array.
\param[in] da   Scalar for the first array.
\param[in] dx   First array.
\param[in] incx Increment of the values in the given first array.
\param[in] dy   Second array.
\param[in] incy Increment of the values in the given second array.

\return DAXPY constant times a vector plus a vector.
*/
void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
#else
/*!
\fn saxpy_

\brief SAXPY constant times a vector plus a vector.

SAXPY constant times a vector plus a vector. Uses unrolled loops for increments
equal to one.

\sa http://www.netlib.org/lapack/explore-html/d8/daf/saxpy_8f.html

\param[in] n    Length of the input array.
\param[in] sa   Scalar for the first array.
\param[in] sx   First array.
\param[in] incx Increment of the values in the given first array.
\param[in] sy   Second array.
\param[in] incy Increment of the values in the given second array.

\return SAXPY constant times a vector plus a vector.
*/
void saxpy_(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\fn dgemv_

\brief Double-precision GEneral matrix Matrix-Vector multiplier.

Performs one of the following matrix-vector operations:

y := alpha*a*x + beta*y, or
y := alpha*a'*x + beta*y.

\sa http://www.math.utah.edu/software/lapack/lapack-blas/dgemv.html

\param[in]      trans   Is this the transpose of the matrix?
\param[in]      m       The number of rows of the matrix a.  m >= 0.
\param[in]      n       The number of columns of the matrix a. n >= 0.
\param[in]      alpha   The scalar alpha.
\param[in,out]  a       The leading m by n part of the array.
\param[in]      lda     The leading dimension of a. lda >= max(1,m).
\param[in,out]  x       Array of DIMENSION at least:
                        (1 + (n - 1)*abs(incx)), and (1 + (m - 1)*abs(incx)), if
                        trans = 'N' or 'n'.
\param[in]      incx    The increment for the elements of x. incx > 0.
\param[in]      beta    The scalar beta.
\param[in,out]  y       Array of DIMENSION at least:
                        (1 + (m - 1)*abs(incy)), and (1 + (n - 1)*abs(incy)), if
                        trans = 'N' or 'n'.
\param[in]      incy    The increment for the elements of y. incy > 0.
*/
void dgemv_(char *trans,
            int *m,
            int *n,
            double *alpha,
            double *a,
            int *lda,
            double *x,
            int *incx,
            double *beta,
            double *y,
            int *incy);
#else
/*!
\fn fgemv_

\brief Single-precision GEneral matrix Matrix-Vector multiplier.

Performs one of the following matrix-vector operations:

y := alpha*a*x + beta*y, or
y := alpha*a'*x + beta*y.

\sa http://www.math.utah.edu/software/lapack/lapack-blas/sgemv.html

\param[in]      trans   Is this the transpose of the matrix?
\param[in]      m       The number of rows of the matrix a.  m >= 0.
\param[in]      n       The number of columns of the matrix a. n >= 0.
\param[in]      alpha   The scalar alpha.
\param[in,out]  a       The leading m by n part of the array.
\param[in]      lda     The leading dimension of a. lda >= max(1,m).
\param[in,out]  x       Array of DIMENSION at least:
                        (1 + (n - 1)*abs(incx)), and (1 + (m - 1)*abs(incx)), if
                        trans = 'N' or 'n'.
\param[in]      incx    The increment for the elements of x. x > 0.
\param[in]      beta    The scalar beta.
\param[in,out]  y       Array of DIMENSION at least:
                        (1 + (m - 1)*abs(incy)), and (1 + (n - 1)*abs(incy)), if
                        trans = 'N' or 'n'.
\param[in]      incy    The increment for the elements of y. incy > 0.
*/
void sgemv_(char *trans,
            int *m,
            int *n,
            float *alpha,
            float *a,
            int *lda,
            float *x,
            int *incx,
            float *beta,
            float *y,
            int *incy);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\fn dgemv_

\brief Double-precision GEneral matrices Matrix-Matrix multiplier.

Performs:

C := alpha*op(a)*op(b) + beta*c

\sa http://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html

\param[in]      transa  Is this the transpose of the matrix a?
\param[in]      transb  Is this the transpose of the matrix b?
\param[in]      m       The number of rows of the matrices a and c.  m >= 0.
\param[in]      n       The number of cols of the matrices b and c.  n >= 0.
\param[in]      k       The number of cols of a and rows of c.  k >= 0.
\param[in]      alpha   Real value used to scale the product of matrices.
\param[in,out]  a       Matrix a.
\param[in]      lda     The leading dimension of a.
\param[in,out]  b       Matrix b.
\param[in]      ldb     The leading dimension of b.
\param[in]      beta    Real value used to scale matrix c.
\param[in,out]  c       Matrix c.
\param[in]      ldc     The leading dimension of c.
*/
void dgemm_(char *transa,
            char* transb,
            int *m,
            int *n,
            int *k,
            double *alpha,
            double *a,
            int *lda,
            double *b,
            int *ldb,
            double *beta,
            double *c,
            int *ldc);
}
#else
/*!
\fn sgemv_

\brief Single-precision GEneral matrices Matrix-Matrix multiplier.

Performs:

C := alpha*op(a)*op(b) + beta*c

\sa http://www.math.utah.edu/software/lapack/lapack-blas/sgemm.html

\param[in]      transa  Is this the transpose of the matrix a?
\param[in]      transb  Is this the transpose of the matrix b?
\param[in]      m       The number of rows of the matrices a and c.  m >= 0.
\param[in]      n       The number of cols of the matrices b and c.  n >= 0.
\param[in]      k       The number of cols of a and rows of c.  k >= 0.
\param[in]      alpha   Real value used to scale the product of matrices.
\param[in,out]  a       Matrix a.
\param[in]      lda     The leading dimension of a.
\param[in,out]  b       Matrix b.
\param[in]      ldb     The leading dimension of b.
\param[in]      beta    Real value used to scale matrix c.
\param[in,out]  c       Matrix c.
\param[in]      ldc     The leading dimension of c.
*/
void sgemm_(char *transa,
            char* transb,
            int *m,
            int *n,
            int *k,
            double *alpha,
            double *a,
            int *lda,
            double *b,aamm
            int *ldb,
            double *beta,
            double *c,
            int *ldc);
}
#endif
}

mtk::Real mtk::BLASAdapter::RealNRM2(Real *in, int &in_length) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(in_length <= 0, __FILE__, __LINE__, __func__);
  #endif

  int incx{1};  // Increment for the elements of xx. ix >= 0.

  #ifdef MTK_PRECISION_DOUBLE
  return dnrm2_(&in_length, in, &incx);
  #else
  return snrm2_(&in_length, in, &incx);
  #endif
}

void mtk::BLASAdapter::RealAXPY(mtk::Real alpha,
                                     mtk::Real *xx,
                                     mtk::Real *yy,
                                     int &in_length) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(xx == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(yy == nullptr, __FILE__, __LINE__, __func__);
  #endif

  int incx{1};  // Increment for the elements of xx. ix >= 0.

  #ifdef MTK_PRECISION_DOUBLE
  daxpy_(&in_length, &alpha, xx, &incx, yy, &incx);
  #else
  saxpy_(&in_length, &alpha, xx, &incx, yy, &incx);
  #endif
}

mtk::Real mtk::BLASAdapter::RelNorm2Error(mtk::Real *computed,
                                          mtk::Real *known,
                                          int length) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(computed == nullptr, __FILE__, __LINE__, __func__);
  mtk::Tools::Prevent(known == nullptr, __FILE__, __LINE__, __func__);
  #endif

  mtk::Real norm_2_computed{mtk::BLASAdapter::RealNRM2(known, length)};

  mtk::Real alpha{-mtk::kOne};

  mtk::BLASAdapter::RealAXPY(alpha, known, computed, length);

  mtk::Real norm_2_difference{mtk::BLASAdapter::RealNRM2(computed, length)};

  return norm_2_difference/norm_2_computed;
}

void mtk::BLASAdapter::RealDenseMV(mtk::Real &alpha,
                                   mtk::DenseMatrix &aa,
                                   mtk::Real *xx,
                                   mtk::Real &beta,
                                   mtk::Real *yy) {

  // Make sure input matrices are row-major ordered.

  if (aa.matrix_properties().ordering() == mtk::COL_MAJOR) {
    aa.OrderRowMajor();
  }

  char transa{'T'}; // State that now, the input WILL be in row-major ordering.

  int mm{aa.num_rows()};                  // Rows of aa.
  int nn{aa.num_cols()};                  // Columns of aa.
  int lda{(aa.matrix_properties()).ld()}; // Leading dimension.
  int incx{1};                            // Increment of values in x.
  int incy{1};                            // Increment of values in y.

  std::swap(mm,nn);
  #ifdef MTK_PRECISION_DOUBLE
  dgemv_(&transa, &mm, &nn, &alpha, aa.data(), &lda,
         xx, &incx, &beta, yy, &incy);
  #else
  sgemv_(&transa, &mm, &nn, &alpha, aa.data(), &lda,
        xx, &incx, &beta, yy, &incy);
  #endif
  std::swap(mm,nn);
}

mtk::DenseMatrix mtk::BLASAdapter::RealDenseMM(mtk::DenseMatrix &aa,
                                               mtk::DenseMatrix &bb) {

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(aa.num_cols() != bb.num_rows(),
                      __FILE__, __LINE__, __func__);
  #endif

  // Make sure input matrices are row-major ordered.

  if (aa.matrix_properties().ordering() == mtk::COL_MAJOR) {
    aa.OrderRowMajor();
  }
  if (bb.matrix_properties().ordering() == mtk::COL_MAJOR) {
    bb.OrderRowMajor();
  }

  char ta{'T'}; // State that input matrix aa is in row-wise ordering.
  char tb{'T'}; // State that input matrix bb is in row-wise ordering.

  int mm{aa.num_rows()};  // Rows of aa and rows of cc.
  int nn{bb.num_cols()};  // Cols of bb and cols of cc.
  int kk{aa.num_cols()};  // Cols of aa and rows of bb.

  int cc_num_rows{mm};  // Rows of cc.
  int cc_num_cols{nn};  // Columns of cc.

  int lda{std::max(1,kk)};  // Leading dimension of the aa matrix.
  int ldb{std::max(1,nn)};  // Leading dimension of the bb matrix.
  int ldc{std::max(1,mm)};  // Leading dimension of the cc matrix.

  mtk::Real alpha{1.0}; // First scalar coefficient.
  mtk::Real beta{0.0};  // Second scalar coefficient.

  mtk::DenseMatrix cc_col_maj_ord(cc_num_rows,cc_num_cols); // Output matrix.

  cc_col_maj_ord.SetOrdering(mtk::COL_MAJOR);

  #ifdef MTK_PRECISION_DOUBLE
  dgemm_(&ta, &tb, &mm, &nn, &kk, &alpha, aa.data(), &lda,
         bb.data(), &ldb, &beta, cc_col_maj_ord.data(), &ldc);
  #else
  sgemm_(&ta, &tb, &mm, &nn, &kk, &alpha, aa.data(), &lda,
         bb.data(), &ldb, &beta, cc_col_maj_ord.data(), &ldc);
  #endif

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "cc_col_maj_ord =" << std::endl;
  std::cout << cc_col_maj_ord << std::endl;
  #endif

  cc_col_maj_ord.OrderRowMajor();

  return cc_col_maj_ord;
}
