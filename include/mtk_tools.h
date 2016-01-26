/*!
\file mtk_tools.h

\brief Tool manager class.

Definition of a class providing basic tools to ensure execution correctness,
and to assists with unitary testing.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\note Performance Tip 8.1. If they do not need to be modified by the called
function, pass large objects using pointers to constant data or references to
constant data, to obtain the performance benefits of pass-by-reference.
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

#ifndef MTK_INCLUDE_TOOLS_H_
#define MTK_INCLUDE_TOOLS_H_

#include <ctime>

#include "mtk_roots.h"

namespace mtk {

/*!
\class Tools

\ingroup c03-execution_tools

\brief Tool manager class.

Basic tools to ensure execution correctness, and to assists with unitary
testing.
*/
class Tools {
 public:
  /*!
  \brief Enforces preconditions by preventing their complements from occur.

  \sa http://stackoverflow.com/questions/8884335/print-the-file-name-line-number-and-function-name-of-a-calling-function-c-pro

  \param[in] complement Complement of desired pre-condition.
  \param[in] fname Name of the file being checked.
  \param[in] lineno Number of the line where the check is executed.
  \param[in] fxname Name of the module containing the check.
  */
  static void Prevent(const bool complement,
                      const char *const fname,
                      int lineno,
                      const char *const fxname) noexcept;

  /*!
  \brief Begins the execution of a unit test. Starts a timer.

  \param[in] nn Number of the test.
  */
  static void BeginUnitTestNo(const int &nn) noexcept;

  /*!
  \brief Ends the execution of a unit test. Stops and reports wall-clock time.

  \param[in] nn Number of the test.
  */
  static void EndUnitTestNo(const int &nn) noexcept;

  /*!
  \brief Asserts if the condition required to pass the unit test occurs.

  \param[in] condition Condition to be asserted.
  */
  static void Assert(const bool &condition) noexcept;

 private:
  static int test_number_;  ///< Current test being executed.

  static Real duration_;    ///< Duration of the current test.

  static clock_t begin_time_;   ///< Elapsed time on current test.
};
}
#endif  // End of: MTK_INCLUDE_TOOLS_H_
