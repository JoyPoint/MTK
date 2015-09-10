/*!
\file mtk_lapack_adapter.cc

\brief Adapter class for the LAPACK API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the LAPACK.

The **LAPACK** is written in Fortran 90 and provides routines for solving
systems of simultaneous linear equations, least-squares solutions of linear
systems of equations, eigenvalue problems, and singular value problems.

\sa http://www.netlib.org/lapack/

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

#include <cstring>

#include <iostream>
#include <iomanip>

#include <algorithm>

#include "mtk_tools.h"
#include "mtk_dense_matrix.h"
#include "mtk_lapack_adapter.h"

namespace mtk {

extern "C" {

#ifdef MTK_PRECISION_DOUBLE
/*!
\fn dgesv_

\brief Double-precision GEneral matrix Simple driVer.

Solves a multiple-rhs dense linear system through LU factorization with partial
pivoting. Matrices have to be in column-major order.

\sa http://www.netlib.org/lapack/explore-html/d8/d72/dgesv_8f.html

\param[in]      n     The order of the matrix A.
\param[in]      nrhs  The number of rows of the matrix a.  m >= 0.
\param[in,out]  a     On entry, the n-by-n matrix a.
\param[in]      lda   The leading dimension of a. lda >= max(1,m).
\param[in,out]  ipiv  The pivots that define the permutation matrix P.
\param[in,out]  b     On entry, matrix b of right-hand side vectors.
\param[in]      ldb   The leading dimension of b. ldb >= max(1,m,n).
\param[in,out]  info  If info = 0, then successful exit.
*/
void dgesv_(int* n,
            int* nrhs,
            Real* a,
            int* lda,
            int* ipiv,
            Real* b,
            int* ldb,
            int* info);
#else
/*!
\fn sgesv_

\brief Single-precision GEneral matrix Simple driVer.

Solves a multiple-rhs dense linear system through LU factorization with partial
pivoting. Matrices have to be in column-major order.

\sa http://www.netlib.org/lapack/explore-html/d7/de8/sgesv_8f.html

\param[in]      n     The order of the matrix A.
\param[in]      nrhs  The number of rows of the matrix a.  m >= 0.
\param[in,out]  a     On entry, the n-by-n matrix a.
\param[in]      lda   The leading dimension of a. lda >= max(1,m).
\param[in,out]  ipiv  The pivots that define the permutation matrix P.
\param[in,out]  b     On entry, matrix b of right-hand side vectors.
\param[in]      ldb   The leading dimension of b. ldb >= max(1,m,n).
\param[in,out]  info  If info = 0, then successful exit.
*/
void sgesv_(int* n,
            int* nrhs,
            Real* a,
            int* lda,
            int* ipiv,
            Real* b,
            int* ldb,
            int* info);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\brief Double-precision GEneral matrix Least Squares solver.

DGELS solves overdetermined or underdetermined real linear systems involving an
M-by-N matrix A, or its transpose, using a QR or LQ factorization of A.  It is
assumed that A has full rank.

The following options are provided:

1. If TRANS = 'N' and m >= n:  find the least squares solution of an
overdetermined system, i.e., solve the least squares problem

                minimize || B - A*X ||.

2. If TRANS = 'N' and m < n:  find the minimum norm solution of an
underdetermined system A * X = B.

3. If TRANS = 'T' and m >= n:  find the minimum norm solution of an
undetermined system A**T * X = B.

4. If TRANS = 'T' and m < n:  find the least squares solution of an
overdetermined system, i.e., solve the least squares problem

                minimize || B - A**T * X ||.

Several right hand side vectors b and solution vectors x can be handled in a
single call; they are stored as the columns of the M-by-NRHS right hand side
matrix B and the N-by-NRHS solution matrix X.

\sa http://www.math.utah.edu/software/lapack/lapack-d/dgels.html

\param[in]      trans Am I giving the transpose of the matrix?
\param[in]      m     The number of rows of the matrix a.  m >= 0.
\param[in]      n     The number of columns of the matrix a.  n >= 0.
\param[in]      nrhs  The number of right-hand sides.
\param[in,out]  a     On entry, the m-by-n matrix a.
\param[in]      lda   The leading dimension of a. lda >= max(1,m).
\param[in,out]  b     On entry, matrix b of right-hand side vectors.
\param[in]      ldb   The leading dimension of b. ldb >= max(1,m,n).
\param[in,out]  work  On exit, if info = 0, work(1) is optimal lwork.
\param[in,out]  lwork The dimension of the array work.
\param[in,out]  info  If info = 0, then successful exit.
*/
void dgels_(char* trans,
            int* m,
            int* n,
            int* nrhs,
            Real* a,
            int* lda,
            Real* b,
            int* ldb,
            Real* work,
            int* lwork,
            int* info);
#else
/*!
\brief Single-precision GEneral matrix Least Squares solver.

SGELS solves overdetermined or underdetermined real linear systems involving an
M-by-N matrix A, or its transpose, using a QR or LQ factorization of A.  It is
assumed that A has full rank.

The following options are provided:

1. If TRANS = 'N' and m >= n:  find the least squares solution of an
overdetermined system, i.e., solve the least squares problem

                minimize || B - A*X ||.

2. If TRANS = 'N' and m < n:  find the minimum norm solution of an
underdetermined system A * X = B.

3. If TRANS = 'T' and m >= n:  find the minimum norm solution of an
undetermined system A**T * X = B.

4. If TRANS = 'T' and m < n:  find the least squares solution of an
overdetermined system, i.e., solve the least squares problem

                minimize || B - A**T * X ||.

Several right hand side vectors b and solution vectors x can be handled in a
single call; they are stored as the columns of the M-by-NRHS right hand side
matrix B and the N-by-NRHS solution matrix X.

\sa http://www.math.utah.edu/software/lapack/lapack-s/sgels.html

\param[in]      trans Am I giving the transpose of the matrix?
\param[in]      m     The number of rows of the matrix a.  m >= 0.
\param[in]      n     The number of columns of the matrix a.  n >= 0.
\param[in]      nrhs  The number of right-hand sides.
\param[in,out]  a     On entry, the m-by-n matrix a.
\param[in]      lda   The leading dimension of a. lda >= max(1,m).
\param[in,out]  b     On entry, matrix b of right-hand side vectors.
\param[in]      ldb   The leading dimension of b. ldb >= max(1,m,n).
\param[in,out]  work  On exit, if info = 0, work(1) is optimal lwork.
\param[in,out]  lwork The dimension of the array work.
\param[in,out]  info  If info = 0, then successful exit.
*/
void sgels_(char* trans,
            int* m,
            int* n,
            int* nrhs,
            Real* a,
            int* lda,
            Real* b,
            int* ldb,
            Real* work,
            int* lwork,
            int* info);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\brief Double-precision GEneral matrix QR Factorization.

Double-Precision Orthogonal Make Q from QR: dormqr_ overwrites the general real
M-by-N matrix C with (Table 1):

                SIDE = 'L'     SIDE = 'R'
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q**T * C       C * Q**T

where Q is a real orthogonal matrix defined as the product of k elementary
reflectors

      Q = H(1) H(2) . . . H(k)

as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N if SIDE =
'R'.

\sa http://www.netlib.org/lapack/lug/node69.html

\param[in]      m     The number of columns of the matrix a.  n >= 0.
\param[in]      n     The number of columns of the matrix a.  n >= 0.
\param[in,out]  a     On entry, the n-by-n matrix a.
\param[in]      lda   Leading dimension matrix.  LDA >= max(1,M).
\param[in,out]  tau   Scalars from elementary reflectors. min(M,N).
\param[in,out]  work  Workspace. info = 0, work(1) is optimal lwork.
\param[in]      lwork The dimension of work. lwork >= max(1,n).
\param[in]      info  info = 0: successful exit.
*/
void dgeqrf_(int *m,
             int *n,
             Real *a,
             int *lda,
             Real *tau,
             Real *work,
             int *lwork,
             int *info);
#else
/*!
\brief Single-precision GEneral matrix QR Factorization.

Single-Precision Orthogonal Make Q from QR: dormqr_ overwrites the general real
M-by-N matrix C with (Table 1):

                SIDE = 'L'     SIDE = 'R'
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q**T * C       C * Q**T

where Q is a real orthogonal matrix defined as the product of k elementary
reflectors

      Q = H(1) H(2) . . . H(k)

as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N if SIDE =
'R'.

\sa http://www.netlib.org/lapack/explore-html/df/d97/sgeqrf_8f.html

\param[in]      m     The number of columns of the matrix a.  n >= 0.
\param[in]      n     The number of columns of the matrix a.  n >= 0.
\param[in,out]  a     On entry, the n-by-n matrix a.
\param[in]      lda   Leading dimension matrix.  LDA >= max(1,M).
\param[in,out]  tau   Scalars from elementary reflectors. min(M,N).
\param[in,out]  work  Workspace. info = 0, work(1) is optimal lwork.
\param[in]      lwork The dimension of work. lwork >= max(1,n).
\param[in]      info  info = 0: successful exit.
*/
void sgeqrf_(int *m,
             int *n,
             Real *a,
             int *lda,
             Real *tau,
             Real *work,
             int *lwork,
             int *info);
#endif

#ifdef MTK_PRECISION_DOUBLE
/*!
\brief Double-precision Orthogonal Matrix from QR factorization.

Double-Precision Orthogonal Make Q from QR: dormqr_ overwrites the general real
M-by-N matrix C with (Table 1):

                SIDE = 'L'     SIDE = 'R'
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q**T * C       C * Q**T

where Q is a real orthogonal matrix defined as the product of k elementary
reflectors

      Q = H(1) H(2) . . . H(k)

as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N if SIDE =
'R'.

\sa http://www.netlib.no/netlib/lapack/Real/dormqr.f

\param[in]      side  See Table 1 above.
\param[in]      trans See Table 1 above.
\param[in]      m     Number of rows of the C matrix.
\param[in]      n     Number of columns of the C matrix.
\param[in]      k     Number of reflectors.
\param[in,out]  a     The matrix containing the reflectors.
\param[in]      lda   The dimension of work. lwork >= max(1,n).
\param[in]      tau   Scalar factors of the elementary reflectors.
\param[in]      c     Output matrix.
\param[in]      ldc   Leading dimension of the output matrix.
\param[in,out]  work  Workspace. info = 0, work(1) optimal lwork.
\param[in]      lwork The dimension of work.
\param[in,out]  info  info = 0: successful exit.
*/
void dormqr_(char *side,
             char *trans,
             int *m,
             int *n,
             int *k,
             Real *a,
             int *lda,
             Real *tau,
             Real *c,
             int *ldc,
             Real *work,
             int *lwork,
             int *info);
#else
/*!
\brief Single-precision Orthogonal Matrix from QR factorization.

Single-Precision Orthogonal Make Q from QR: sormqr_ overwrites the general real
M-by-N matrix C with (Table 1):

                SIDE = 'L'     SIDE = 'R'
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q**T * C       C * Q**T

where Q is a real orthogonal matrix defined as the product of k elementary
reflectors

      Q = H(1) H(2) . . . H(k)

as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N if SIDE =
'R'.

\sa http://www.netlib.org/lapack/explore-html/d0/d98/sormqr_8f_source.html

\param[in]      side  See Table 1 above.
\param[in]      trans See Table 1 above.
\param[in]      m     Number of rows of the C matrix.
\param[in]      n     Number of columns of the C matrix.
\param[in]      k     Number of reflectors.
\param[in,out]  a     The matrix containing the reflectors.
\param[in]      lda   The dimension of work. lwork >= max(1,n).
\param[in]      tau   Scalar factors of the elementary reflectors.
\param[in]      c     Output matrix.
\param[in]      ldc   Leading dimension of the output matrix.
\param[in,out]  work  Workspace. info = 0, work(1) optimal lwork.
\param[in]      lwork The dimension of work.
\param[in,out]  info  info = 0: successful exit.
*/
void sormqr_(char *side,
             char *trans,
             int *m,
             int *n,
             int *k,
             Real *a,
             int *lda,
             Real *tau,
             Real *c,
             int *ldc,
             Real *work,
             int *lwork,
             int *info);
#endif
}
}

int mtk::LAPACKAdapter::SolveDenseSystem(mtk::DenseMatrix &mm,
                                         mtk::Real *rhs) {


  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(rhs == nullptr, __FILE__, __LINE__, __func__);
  #endif

  int *ipiv{};                // Array for pivoting information.
  int nrhs{1};                // Number of right-hand sides.
  int info{};                 // Status of the solution.
  int mm_rank{mm.num_rows()}; // Rank of the matrix.

  try {
    ipiv = new int[mm_rank];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(ipiv, 0, sizeof(ipiv[0])*mm_rank);

  int ldbb = mm_rank;
  int mm_ld = mm_rank;

  #ifdef MTK_PRECISION_DOUBLE
  dgesv_(&mm_rank, &nrhs, mm.data(), &mm_ld, ipiv, rhs, &ldbb, &info);
  #else
  fgesv_(&mm_rank, &nrhs, mm.data(), &mm_ld, ipiv, rhs, &ldbb, &info);
  #endif

  delete [] ipiv;

  return info;
}

int mtk::LAPACKAdapter::SolveDenseSystem(mtk::DenseMatrix &mm,
                                         mtk::DenseMatrix &bb) {

  int nrhs{bb.num_rows()};  // Number of right-hand sides.

  #if MTK_DEBUG_LEVEL > 0
  mtk::Tools::Prevent(nrhs <= 0, __FILE__, __LINE__, __func__);
  #endif

  int *ipiv{};                // Array for pivoting information.
  int info{};                 // Status of the solution.
  int mm_rank{mm.num_rows()}; // Rank of the matrix.

  try {
    ipiv = new int[mm_rank];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(ipiv, 0, sizeof(ipiv[0])*mm_rank);

  int ldbb = mm_rank;
  int mm_ld = mm_rank;

  #ifdef MTK_PRECISION_DOUBLE
  dgesv_(&mm_rank, &nrhs, mm.data(), &mm_ld, ipiv, bb.data(), &ldbb, &info);
  #else
  fgesv_(&mm_rank, &nrhs, mm.data(), &mm_ld, ipiv, bb.data(), &ldbb, &info);
  #endif

  delete [] ipiv;

  // After output, the data in the matrix will be column-major ordered.

  bb.SetOrdering(mtk::COL_MAJOR);

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "bb_col_maj_ord =" << std::endl;
  std::cout << bb << std::endl;
  #endif

  bb.OrderRowMajor();

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "bb_row_maj_ord =" << std::endl;
  std::cout << bb << std::endl;
  #endif

  return info;
}

mtk::DenseMatrix mtk::LAPACKAdapter::QRFactorDenseMatrix(mtk::DenseMatrix &aa) {

  mtk::Real *work{}; // Working array.
  mtk::Real *tau{};  // Array for the Householder scalars.

  // Prepare to factorize: allocate and inquire for the value of lwork.
  try {
    work = new mtk::Real[1];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(work, mtk::kZero, sizeof(aa.data()[0])*1);

  int lwork{-1};
  int info{};

  int aa_num_cols = aa.num_cols();
  int aaT_num_rows = aa.num_cols();
  int aaT_num_cols = aa.num_rows();

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Input matrix BEFORE QR factorization:" << std::endl;
  std::cout << aa << std::endl;
  #endif

  #ifdef MTK_PRECISION_DOUBLE
  dgeqrf_(&aaT_num_rows, &aaT_num_cols, aa.data(), &aaT_num_rows,
          tau,
          work, &lwork, &info);
  #else
  fgeqrf_(&aaT_num_rows, &aaT_num_cols, aa.data(), &aaT_num_rows,
          tau,
          work, &lwork, &info);
  #endif

  #if MTK_DEBUG_LEVEL > 0
  if (info == 0) {
    lwork = (int) work[0];
  } else {
    std::cerr << "Could not get value for lwork on line " << __LINE__ - 5 <<
      std::endl;
    std::cerr << "Exiting..." << std::endl;
  }
  #endif

  #if MTK_DEBUG_LEVEL>0
  std::cout << "lwork = " << std::endl << std::setw(12) << lwork << std::endl
    << std::endl;
  #endif

  delete [] work;
  work = nullptr;

  // Once we know lwork, we can actually invoke the factorization:
  try {
    work = new mtk::Real [lwork];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(work, mtk::kZero, sizeof(work[0])*lwork);

  int ltau = std::min(aaT_num_rows,aaT_num_cols);

  try {
    tau = new mtk::Real [ltau];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(tau, mtk::kZero, sizeof(0.0)*ltau);

  #ifdef MTK_PRECISION_DOUBLE
  dgeqrf_(&aaT_num_rows, &aaT_num_cols, aa.data(), &aaT_num_rows,
          tau, work, &lwork, &info);
  #else
  fgeqrf_(&aaT_num_rows, &aaT_num_cols, aa.data(), &aaT_num_rows,
          tau, work, &lwork, &info);
  #endif

  if (!info) {
    #if MTK_DEBUG_LEVEL > 0
    std::cout << "QR factorization completed!" << std::endl << std::endl;
    #endif
  } else {
    std::cerr << "Error solving system! info = " << info << std::endl;
    std::cerr << "Exiting..." << std::endl;
  }

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Input matrix AFTER QR factorization:" << std::endl;
  std::cout << aa << std::endl;
  #endif

  // We now generate the real matrix Q with orthonormal columns. This has to
  // be done separately since the actual output of dgeqrf_ (AA_) represents
  // the orthogonal matrix Q as a product of min(aa_num_rows,aa_num_cols)
  // elementary Householder reflectors. Notice that we must re-inquire the new
  // value for lwork that is used.

  bool padded{false};

  bool transpose{false};

  mtk::DenseMatrix QQ_(aa.num_cols(),padded,transpose);

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Initialized QQ_T: " << std::endl;
  std::cout << QQ_ << std::endl;
  #endif

  // Assemble the QQ_ matrix:
  lwork = -1;

  delete[] work;
  work = nullptr;

  try {
    work = new mtk::Real[1];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() <<
      std::endl;
  }
  memset(work, mtk::kZero, sizeof(work[0])*1);

  char side_{'L'};
  char trans_{'N'};

  int aux = QQ_.num_rows();

  #ifdef MTK_PRECISION_DOUBLE
  dormqr_(&side_, &trans_,
          &aa_num_cols, &aa_num_cols, &ltau, aa.data(), &aaT_num_rows, tau,
          QQ_.data(), &aux, work, &lwork, &info);
  #else
  formqr_(&side_, &trans_,
          &aa_num_cols, &aa_num_cols, &ltau, aa.data(), &aaT_num_rows, tau,
          QQ_.data(), &aux, work, &lwork, &info);
  #endif

  #if MTK_DEBUG_LEVEL > 0
  if (info == 0) {
    lwork = (int) work[0];
  } else {
    std::cerr << "Could not get lwork on line " << __LINE__ - 5 << std::endl;
    std::cerr << "Exiting..." << std::endl;
  }
  #endif

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "lwork = " << std::endl << std::setw(12) << lwork <<
    std::endl << std::endl;
  #endif

  delete[] work;
  work = nullptr;

  try {
    work = new mtk::Real[lwork];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(work, mtk::kZero, sizeof(work[0])*lwork);

  #ifdef MTK_PRECISION_DOUBLE
  dormqr_(&side_, &trans_,
          &aa_num_cols, &aa_num_cols, &ltau, aa.data(), &aaT_num_rows, tau,
          QQ_.data(), &aux, work, &lwork, &info);
  #else
  formqr_(&side_, &trans_,
          &aa_num_cols, &aa_num_cols, &ltau, aa.data(), &aaT_num_rows, tau,
          QQ_.data(), &aux, work, &lwork, &info);
  #endif

  if (!info) {
    #if MTK_DEBUG_LEVEL>0
    std::cout << "Q matrix successfully assembled!" << std::endl << std::endl;
    #endif
  } else {
    std::cerr << "Something went wrong solving system! info = " << info <<
      std::endl;
    std::cerr << "Exiting..." << std::endl;
  }

  delete[] work;
  work = nullptr;

  delete[] tau;
  tau = nullptr;

  return QQ_;
}

int mtk::LAPACKAdapter::SolveRectangularDenseSystem(const mtk::DenseMatrix &aa,
                                                   mtk::Real *ob_,
                                                   int ob_ld_) {

  // We first invoke the solver to query for the value of lwork. For this,
  // we must at least allocate enough space to allow access to WORK(1), or
  // work[0]:

  // If LWORK = -1, then a workspace query is assumed; the routine only
  // calculates the optimal size of the WORK array, returns this value as
  // the first entry of the WORK array, and no error message related to
  // LWORK is issued by XERBLA.

  mtk::Real *work{}; // Work array.

  try {
    work = new mtk::Real[1];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 << std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(work, mtk::kZero, sizeof(work[0])*1);

  char trans_{'N'};
  int nrhs_{1};
  int info{0};
  int lwork{-1};

  int AA_num_rows_  = aa.num_cols();
  int AA_num_cols_  = aa.num_rows();
  int AA_ld_ = std::max(1,aa.num_cols());

  #ifdef MTK_PRECISION_DOUBLE
  dgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, aa.data(), &AA_ld_,
         ob_, &ob_ld_,
         work, &lwork, &info);
  #else
  sgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, aa.data(), &AA_ld_,
         ob_, &ob_ld_,
         work, &lwork, &info);
  #endif

  if (info == 0) {
    lwork = (int) work[0];
  } else {
    std::cerr << "Could not get value for lwork on line " << __LINE__ - 2 <<
      std::endl;
    std::cerr << "Exiting..." << std::endl;
    return info;
  }

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "lwork = " << std::endl << std::setw(12)<< lwork <<
    std::endl << std::endl;
  #endif

  // We then use lwork's new value to create the work array:
  delete[] work;
  work = nullptr;

  try {
    work = new mtk::Real[lwork];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 << std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(work, 0.0, sizeof(work[0])*lwork);

  // We now invoke the solver again:
  #ifdef MTK_PRECISION_DOUBLE
  dgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, aa.data(), &AA_ld_,
         ob_, &ob_ld_,
         work, &lwork, &info);
  #else
  sgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, aa.data(), &AA_ld_,
         ob_, &ob_ld_,
         work, &lwork, &info);
  #endif

  delete [] work;
  work = nullptr;

  return info;
}
