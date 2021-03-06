/*!
\file 1d_positive_weights.cc

\brief The CBS algorithm computes positive-definite weights, for 1D operators.

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

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>

#include <vector>

#include "mtk.h"

int main () {

  std::cout << "Example: Positive-Definite Weights for 1D Mimetic"
    "Operators." << std::endl;

  /// 1. Create all critical-order divergence operators.

  mtk::Grad1D grad10;

  bool assertion = grad10.ConstructGrad1D(10);
  if (!assertion) {
    std::cerr << "Mimetic grad (10th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::Grad1D grad12;

  assertion = grad12.ConstructGrad1D(12);
  if (!assertion) {
    std::cerr << "Mimetic grad (12th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::Grad1D grad14;

  assertion = grad14.ConstructGrad1D(14);
  if (!assertion) {
    std::cerr << "Mimetic grad (14th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  /// 2. Create all critical-order divergence operators.

  mtk::Div1D div8;

  assertion = div8.ConstructDiv1D(8);
  if (!assertion) {
    std::cerr << "Mimetic div (8th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::Div1D div10;

  assertion = div10.ConstructDiv1D(10);
  if (!assertion) {
    std::cerr << "Mimetic div (10th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::Div1D div12;

  assertion = div12.ConstructDiv1D(12);
  if (!assertion) {
    std::cerr << "Mimetic div (12th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::Div1D div14;

  assertion = div14.ConstructDiv1D(14);
  if (!assertion) {
    std::cerr << "Mimetic div (14th order) could not be built." << std::endl;
    return EXIT_FAILURE;
  }
}
