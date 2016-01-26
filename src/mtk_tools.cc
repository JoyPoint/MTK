/*!
\file mtk_tools.cc

\brief Tool manager class.

Implementation of a class providing basic tools to ensure execution correctness,
and to assists with unitary testing.

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

#include <iostream>

#include "mtk_roots.h"
#include "mtk_tools.h"

void mtk::Tools::Prevent(const bool condition,
                         const char *const fname,
                         int lineno,
                         const char *const fxname) noexcept {

  /// \todo Check if this is the best way of stalling execution.
  if (lineno < 1) {
    std::cerr << __FILE__ << ": " << "Incorrect parameter at line " <<
    __LINE__ - 2 << " (" << __func__ << ")" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (condition) {
    std::cerr << fname << ": " << "Incorrect parameter at line " <<
    lineno << " (" << fxname << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}

int mtk::Tools::test_number_{};   // Current test being executed.

mtk::Real mtk::Tools::duration_{};    // Duration of the current test.

clock_t mtk::Tools::begin_time_{};  // Elapsed time on current test.

void mtk::Tools::BeginUnitTestNo(const int &nn) noexcept {

  #if MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(nn <= 0, __FILE__, __LINE__, __func__);
  #endif

  test_number_ = nn;

  std::cout << "Beginning test " << nn << "." << std::endl;
  begin_time_ = clock();
}

void mtk::Tools::EndUnitTestNo(const int &nn) noexcept {

  #if MTK_PERFORM_PREVENTIONS
  mtk::Tools::Prevent(test_number_ != nn, __FILE__, __LINE__, __func__);
  #endif

  duration_ = mtk::Real(clock() - begin_time_)/CLOCKS_PER_SEC;
}

void mtk::Tools::Assert(const bool &condition) noexcept {

  if (condition) {
    std::cout << "Test " << test_number_ << ": PASSED in " << duration_ <<
      " s." << std::endl;
  } else {
    std::cout << "Test " << test_number_ << ": FAILED in " << duration_ <<
      " s." << std::endl;
  }
}
