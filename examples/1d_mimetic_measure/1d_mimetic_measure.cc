/*!
\file 1d_mimetic_measure.cc

\brief Compute the mimetic measure of different mimetic operators.



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

#if __cplusplus == 201103L

#include <cmath>

#include <iostream>
#include <fstream>

#include <array>

#include "mtk.h"

int main () {

  std::cout << "Example: Computing the mimetic measure of gradient, "
    "divergence, and Laplacian operators." << std::endl;

  const int max_order{14};

  std::ofstream output_dat_file;  // Output file.

  output_dat_file.open("mimetic_measure.tex");

  if (!output_dat_file.is_open()) {
    std::cerr << "Could not open data file." << std::endl;
    return EXIT_FAILURE;
  }

  output_dat_file << "\\begin{tabular}[c]{c:ccc}" << std::endl;
  output_dat_file << "\\toprule" << std::endl;
  output_dat_file << "$k$ & $\\mu(\\breve{\\mathbf{G}}^k_x)$ &"
    "$\\mu(\\breve{\\mathbf{D}}^k_x)$ & "
    "$\\mu(\\breve{\\mathbf{L}}^k_x)$ \\\\" << std::endl;
  output_dat_file << "\\midrule" << std::endl;

  mtk::Grad1D grad;
  mtk::Div1D div;
  mtk::Lap1D lap;

  for (int order = 2; order <= max_order; order += 2) {

    bool go1{grad.ConstructGrad1D(order)};
    bool go2{div.ConstructDiv1D(order)};
    bool go3{lap.ConstructLap1D(order)};

    if (go1 && go2 && go3) {
      output_dat_file << order << " & " << grad.mimetic_measure() << " & " <<
        div.mimetic_measure() << " & " << lap.mimetic_measure() << "\\\\" <<
        std::endl;

    } else {
      std::cerr << "Mimetic operator could not be built." << std::endl;
      return EXIT_FAILURE;
    }
  }

  output_dat_file << "\\bottomrule" << std::endl;
  output_dat_file << "\\end{tabular}" << std::endl;

  output_dat_file.close();
}
#else
#include <iostream>
using std::cout;
using std::endl;
int main () {
  cout << "This code HAS to be compiled with support for C++11." << endl;
  cout << "Exiting..." << endl;
  return EXIT_SUCCESS;
}
#endif
