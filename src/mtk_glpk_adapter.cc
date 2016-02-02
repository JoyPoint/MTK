/*!
\file mtk_glpk_adapter.cc

\brief Adapter class for the GLPK API.

Implementation of a class that contains a collection of static member
functions, that posses direct access to the underlying structure of the
matrices, thus allowing programmers to exploit some of the numerical methods
implemented in the GLPK.

The **GLPK (GNU Linear Programming Kit)** package is intended for solving
large-scale linear programming (LP), mixed integer programming (MIP), and other
related problems. It is a set of routines written in ANSI C and organized in the
form of a callable library.

\sa http://www.gnu.org/software/glpk/

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

#include <cmath>
#include <cstring>

#include <iostream>
#include <iomanip>
#include <limits>

#include "mtk_roots.h"
#include "mtk_blas_adapter.h"
#include "mtk_glpk_adapter.h"

int mtk::GLPKAdapter::num_feasible_solutions_{};

int mtk::GLPKAdapter::num_feasible_solutions() noexcept {

  return num_feasible_solutions_;
}

mtk::Real mtk::GLPKAdapter::SolveSimplexAndCompare(mtk::Real *A,
                                                   int nrows,
                                                   int ncols,
                                                   int kk,
                                                   mtk::Real *hh,
                                                   mtk::Real *qq,
                                                   int robjective,
                                                   mtk::Real mimetic_threshold,
                                                   int copy) noexcept {

  #if MTK_DEBUG_LEVEL > 0
  char mps_file_name[18]; // File name for the MPS files.
  #endif
  char rname[5];          // Row name.
  char cname[5];          // Column name.

  glp_prob *lp; // Linear programming problem.

  int *ia;  // Array for the problem.
  int *ja;  // Array for the problem.

  int problem_size; // Size of the problem.
  int lp_nrows;     // Number of rows.
  int lp_ncols;     // Number of columns.
  int matsize;      // Size of the matrix.
  int glp_index{1}; // Index of the objective function.
  int ii;           // Iterator.
  int jj;           // Iterator.

  mtk::Real *ar;            // Array for the problem.
  mtk::Real *objective;     // Array containing the objective function.
  mtk::Real *rhs;           // Array containing the rhs.
  mtk::Real *err;           // Array of errors.

  mtk::Real x1;             // Norm-2 of the error.

  #if MTK_DEBUG_LEVEL > 0
  mtk::Real obj_value;      // Value of the objective function.
  #endif

  lp_nrows = kk;
  lp_ncols = kk;

  matsize = lp_nrows*lp_ncols;

  /// \warning GLPK indexes in [1,n], so we must get the extra space needed.

  /// 1. Memory allocation.

  problem_size = lp_nrows*lp_ncols + 1;

  try {
    ia = new int[problem_size];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(ia, 0, sizeof(ia[0])*problem_size);

  try {
    ja = new int[problem_size];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(ja, 0, sizeof(ja[0])*problem_size);

  try {
    ar = new mtk::Real[problem_size];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(ar, mtk::kZero, sizeof(ar[0])*problem_size);

  try {
    objective = new mtk::Real[lp_ncols + 1];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(objective, mtk::kZero, sizeof(objective[0])*(lp_ncols + 1));

  try {
    rhs = new mtk::Real[lp_nrows + 1];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(rhs, mtk::kZero, sizeof(rhs[0])*(lp_nrows + 1));

  try {
    err = new mtk::Real[lp_nrows];
  } catch (std::bad_alloc &memory_allocation_exception) {
    std::cerr << "Memory allocation exception on line " << __LINE__ - 3 <<
      std::endl;
    std::cerr << memory_allocation_exception.what() << std::endl;
  }
  memset(err, mtk::kZero, sizeof(err[0])*(lp_nrows));

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "Problem size: " << problem_size << std::endl;
  std::cout << "lp_nrows = " << lp_nrows << std::endl;
  std::cout << "lp_ncols = " << lp_ncols << std::endl;
  std::cout << std::endl;
  #endif

  lp = glp_create_prob();

  glp_set_prob_name (lp, "mtk::GLPKAdapter::Simplex");

  glp_set_obj_dir (lp, GLP_MIN);

  /// 2. Fill the problem.

  glp_add_rows(lp, lp_nrows);

  for (ii = 1; ii <= lp_nrows; ++ii) {
    sprintf(rname, "R%02d",ii);
    glp_set_row_name(lp, ii, rname);
  }

  glp_add_cols(lp, lp_ncols);

  for (ii = 1; ii <= lp_ncols; ++ii) {
    sprintf(cname, "Q%02d",ii);
    glp_set_col_name (lp, ii, cname);
  }

  /// 3. Copy the row to the vector objective.

  #if MTK_DEBUG_LEVEL>0
  std::cout << "Using row " << robjective + 1 << " as objective." << std::endl;
  #endif
  for (jj = 0; jj < kk; ++jj) {
    objective[glp_index] = A[jj + robjective * ncols];
    glp_index++;
  }
  #if MTK_DEBUG_LEVEL >0
  std::cout << std::endl;
  #endif

  /// 4. Forming the RHS.

  glp_index = 1;
  rhs[0] = mtk::kZero;
  for (ii = 0; ii <= lp_nrows; ++ii) {
    if (ii != robjective) {
      rhs[glp_index] = hh[ii];
      glp_set_row_bnds(lp, glp_index, GLP_UP, 0.0, rhs[glp_index]);
      glp_index++;
    }
  }

  #if MTK_DEBUG_LEVEL > 0
  std::cout << "rhs =" << std::endl;
  for (auto ii = 0; ii < lp_nrows; ++ii) {
    std::cout << std::setw(15) << rhs[ii] << std::endl;
  }
  std::cout << std::endl;
  #endif

  /// 5. Setting up the objective function.

  for (ii = 1; ii <= lp_ncols; ++ii) {
    glp_set_obj_coef (lp, ii, objective[ii]);
  }

  /// 6. Setting up constraints.

  for (ii = 1; ii <= lp_ncols; ++ii) {
    glp_set_col_bnds (lp, ii, GLP_LO, mimetic_threshold, 0.0);
  }

  /// 7. Copy the matrix minus the row objective to the glpk problem.

  glp_index = 1;
  for (ii = 0; ii <= kk; ++ii) {
    for (jj = 0; jj < kk; ++jj) {
      if (ii != robjective) {
        ar[glp_index] = A[jj + ii * ncols];
        glp_index++;
      }
    }
  }

  glp_index = 0;

  for (ii = 1; ii < problem_size; ++ii) {
    if (((ii - 1) % lp_ncols) == 0) {
      glp_index++;
    }
    ia[ii] = glp_index;
    ja[ii] = (ii - 1) % lp_ncols + 1;
  }

  glp_load_matrix (lp, matsize, ia, ja, ar);

  #if MTK_DEBUG_LEVEL > 0
  sprintf(mps_file_name, "LP_MPS_row_%02d.mps", robjective);
  glp_write_mps(lp, GLP_MPS_FILE, nullptr, mps_file_name);
  #endif

  /// 8. Solve problem.

  glp_simplex (lp, nullptr);

  // Check status of the solution, determining if this was a feasible solution.

  if (glp_get_status(lp) == GLP_OPT) {

    for(ii = 1; ii <= lp_ncols; ++ii) {
      err[ii - 1] = qq[ii - 1] - glp_get_col_prim(lp,ii);
    }

    #if MTK_DEBUG_LEVEL > 0
    obj_value = glp_get_obj_val (lp);
    std::cout << std::setw(12) << "CBS" << std::setw(12) << "CRS" << std::endl;
    for (ii = 0; ii < lp_ncols; ++ii) {
      std::cout << "q_" << ii + 1 << " = " << std::setw(12) <<
        glp_get_col_prim(lp,ii + 1) << std::setw(12) << qq[ii] << std::endl;
    }
    std::cout << "Objective function value (row " << robjective + 1 << ") = " <<
      obj_value << std::endl;
    #endif

    if (copy) {
      for(ii = 0; ii < lp_ncols; ++ii) {
        qq[ii] = glp_get_col_prim(lp,ii + 1);
      }
      // Preserve the bottom values of qq.
    }

    x1 = mtk::BLASAdapter::RealNRM2(err,lp_ncols);

    num_feasible_solutions_++;

  } else {
    x1 = std::numeric_limits<mtk::Real>::infinity();
  }

  glp_delete_prob (lp);
  glp_free_env ();

  delete [] ia;
  delete [] ja;
  delete [] ar;
  delete [] objective;
  delete [] rhs;
  delete [] err;

  return x1;
}
