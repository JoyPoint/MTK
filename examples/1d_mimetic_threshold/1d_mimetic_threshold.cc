/*!
\file 1d_mimetic_threshold.cc

\brief Perform a sensitivity analysis of the mimetic threshold.



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

#include <cmath>

#include <iostream>
#include <fstream>

#include <array>

#include "mtk.h"

int main () {

  std::cout << "Example: Sensitivity analysis of the mimetic threshold when "
    "constructing 1D gradient and divergence operators." << std::endl;

  std::array<int, 4> orders = {8, 10, 12, 14};

  const int num_samples{10};

  std::vector<mtk::Real> thresholds(num_samples);

  thresholds.at(0) = 1e0;
  for (size_t ii = 1; ii < thresholds.size(); ++ii) {
    thresholds.at(ii) = thresholds.at(ii - 1)/10.0;
  }

  std::array<mtk::Real, num_samples> num_feasible_sols{};

  std::ofstream output_dat_file;  // Output file.

  for (int kk: orders) {

    output_dat_file.open("mimetic_threshold_div_" + std::to_string(kk) +
      ".dat");

    if (!output_dat_file.is_open()) {
      std::cerr << "Could not open data file." << std::endl;
      return EXIT_FAILURE;
    }

    output_dat_file << "# " << 'e' <<  ' ' << 'n' << std::endl;

    for (size_t ii = 0; ii < thresholds.size(); ++ii) {

      mtk::Real epsilon{thresholds.at(ii)};

      mtk::Div1D div;

      if (div.ConstructDiv1D(kk, epsilon)) {
        num_feasible_sols[ii] = div.num_feasible_sols();
      } else {
        std::cerr << "Mimetic divergence could not be built." << std::endl;
        return EXIT_FAILURE;
      }
    }

    for (unsigned int ii = 0; ii < thresholds.size(); ++ii) {
      output_dat_file << thresholds[ii] << ' ' << num_feasible_sols[ii] <<
        std::endl;
    }

    output_dat_file.close();
  }

  for (int kk: orders) {

    output_dat_file.open("mimetic_threshold_grad_" + std::to_string(kk) +
      ".dat");

    if (!output_dat_file.is_open()) {
      std::cerr << "Could not open data file." << std::endl;
      return EXIT_FAILURE;
    }

    output_dat_file << "# " << 'e' <<  ' ' << 'n' << std::endl;

    for (size_t ii = 0; ii < thresholds.size(); ++ii) {

      mtk::Real epsilon{thresholds.at(ii)};

      mtk::Grad1D grad;

      if (grad.ConstructGrad1D(kk, epsilon)) {
        num_feasible_sols[ii] = grad.num_feasible_sols();
      } else {
        std::cerr << "Mimetic gradient could not be built." << std::endl;
        return EXIT_FAILURE;
      }
    }

    for (unsigned int ii = 0; ii < thresholds.size(); ++ii) {
      output_dat_file << thresholds[ii] << ' ' << num_feasible_sols[ii] <<
        std::endl;
    }

    output_dat_file.close();
  }
}
