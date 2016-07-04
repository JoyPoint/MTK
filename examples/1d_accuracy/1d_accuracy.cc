/*!
\file 1d_accuracy.cc

\brief Check the accuracy of mimetic operators.



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

#if __cplusplus == 201103L

#include <cmath>

#include <iostream>
#include <fstream>

#include <array>
#include <vector>

#include "mtk.h"

/// \todo Use Horner's method to evaluate these:
/// https://en.wikipedia.org/wiki/Horner's_method

mtk::Real Polynomial(const mtk::Real &xx, const std::vector<mtk::Real> &pp) {

  mtk::Real sum{};

  mtk::Real kk{pp[0]};

  for (int ii = 0; ii <= kk; ++ii) {
    sum += pow(xx, ii);
  }
  return sum;
}

mtk::Real PolynomialDerivative1(const mtk::Real &xx,
                                const std::vector<mtk::Real> &pp) {

  mtk::Real sum{};

  mtk::Real kk{pp[0]};

  for (int ii = 1; ii <= kk; ++ii) {
    sum += ii*pow(xx, ii - 1);
  }
  return sum;
}

mtk::Real HomogeneousDirichlet(const mtk::Real &tt,
		                           const std::vector<mtk::Real> &pp) {

  return mtk::kOne;
}

int main () {

  std::cout << "Example: Checking the accuracy of the mimetic operators." <<
    std::endl;

  const int max_order{14};

  mtk::Real aa = 0.0;
  mtk::Real bb = 1.0;

  std::array<mtk::Real, 3> epsilons{1.0e-3, 1.0e-6, 1.0e-9};

  std::ofstream output_dat_file;  // Output file.

  /// 1. Perform experiment for the gradient operator.

  output_dat_file.open("accuracy_grad.tex");

  if (!output_dat_file.is_open()) {
    std::cerr << "Could not open data file." << std::endl;
    return EXIT_FAILURE;
  }

  output_dat_file << "\\begin{tabular}[c]{c:ccc}" << std::endl;
  output_dat_file << "\\toprule" << std::endl;
  output_dat_file << "$k$ & \\multicolumn{3}{c}{Relative error} \\\\" <<
    std::endl;
  output_dat_file << " & $\\epsilon = 1\\times 10^{-3}$ &"
    "$\\epsilon = 1\\times 10^{-6}$ & "
    "$\\epsilon = 1\\times 10^{-9}$ \\\\" << std::endl;
  output_dat_file << "\\midrule" << std::endl;

  mtk::Grad1D grad;

  for (int order = 2; order <= max_order; order += 2) {

    output_dat_file << order;

    for (mtk::Real &epsilon: epsilons) {

      if (!grad.ConstructGrad1D(order, epsilon)) {
        std::cerr << "Mimetic gradient could not be built." << std::endl;
        return EXIT_FAILURE;
      }

      int nn{500};

      mtk::UniStgGrid1D gg(aa, bb, nn);

      std::vector<mtk::Real> order_in;

      order_in.push_back((mtk::Real) order);

      gg.BindScalarField(Polynomial, order_in);

      if (!gg.WriteToFile("1d_accuracy_gg.dat", "x", "p(x)")) {
        std::cerr << "Polynomial could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::DenseMatrix gradm{grad.ReturnAsDenseMatrix(gg)};

      if (!gradm.WriteToFile("1d_accuracy_gradm.dat")) {
        std::cerr << "Grad matrix could not be written to disk." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::UniStgGrid1D out(gg.west_bndy_x(),
                          gg.east_bndy_x(),
                          gg.num_cells_x(),
                          mtk::FieldNature::VECTOR);

      out.GenerateDiscreteDomainX();
      out.ReserveDiscreteField();

      mtk::OperatorApplicator::ApplyDenseMatrixGradientOn1DGrid(gradm, gg, out);

      if (!out.WriteToFile("1d_accuracy_out.dat", "x", "~p'(x)")) {
        std::cerr << "Comp Grad Pol could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::UniStgGrid1D gg_prime(aa, bb, nn, mtk::FieldNature::VECTOR);

      gg_prime.BindVectorField(PolynomialDerivative1, order_in);

      if (!gg_prime.WriteToFile("1d_accuracy_gg_prime.dat", "x", "p'(x)")) {
        std::cerr << "Grad Poly could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::Real relative_norm_2_error{};

      relative_norm_2_error =
        mtk::BLASAdapter::RelNorm2Error(out.discrete_field(),
                                        gg_prime.discrete_field(),
                                        gg_prime.num_cells_x());

      output_dat_file <<  " & " << relative_norm_2_error;
    }

    output_dat_file <<  "\\\\" << std::endl;
  }

  output_dat_file << "\\bottomrule" << std::endl;
  output_dat_file << "\\end{tabular}" << std::endl;

  output_dat_file.close();

  /// 2. Perform experiment for the divergence operator.

  output_dat_file.open("accuracy_div.tex");

  if (!output_dat_file.is_open()) {
    std::cerr << "Could not open data file." << std::endl;
    return EXIT_FAILURE;
  }

  output_dat_file << "\\begin{tabular}[c]{c:ccc}" << std::endl;
  output_dat_file << "\\toprule" << std::endl;
  output_dat_file << "$k$ & \\multicolumn{3}{c}{Relative error} \\\\" <<
    std::endl;
  output_dat_file << " & $\\epsilon = 1\\times 10^{-3}$ &"
    "$\\epsilon = 1\\times 10^{-6}$ & "
    "$\\epsilon = 1\\times 10^{-9}$ \\\\" << std::endl;
  output_dat_file << "\\midrule" << std::endl;

  mtk::Div1D div;

  for (int order = 2; order <= max_order; order += 2) {

    output_dat_file << order;

    for (mtk::Real &epsilon: epsilons) {

      if (!div.ConstructDiv1D(order, epsilon)) {
        std::cerr << "Mimetic divergence could not be built." << std::endl;
        return EXIT_FAILURE;
      }

      int nn{500};

      mtk::UniStgGrid1D vv(aa, bb, nn, mtk::FieldNature::VECTOR);

      std::vector<mtk::Real> order_in;

      order_in.push_back((mtk::Real) order);

      vv.BindVectorField(Polynomial, order_in);

      if (!vv.WriteToFile("1d_accuracy_vv.dat", "x", "||p(x)||")) {
        std::cerr << "Polynomial could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::DenseMatrix divm{div.ReturnAsDenseMatrix(vv)};

      mtk::RobinBCDescriptor1D robin_bc_desc_1d;

      robin_bc_desc_1d.PushBackWestCoeff(HomogeneousDirichlet);
      robin_bc_desc_1d.PushBackEastCoeff(HomogeneousDirichlet);

      if (!robin_bc_desc_1d.ImposeOnDivergenceMatrix(div, divm)) {
        std::cerr << "BCs  could not be bound to the matrix." << std::endl;
        return EXIT_FAILURE;
      }

      if (!divm.WriteToFile("1d_accuracy_divm.dat")) {
        std::cerr << "Div matrix could not be written to disk." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::UniStgGrid1D out2(vv.west_bndy_x(),
                             vv.east_bndy_x(),
                             vv.num_cells_x());

      out2.GenerateDiscreteDomainX();
      out2.ReserveDiscreteField();

      mtk::OperatorApplicator::ApplyDenseMatrixDivergenceOn1DGrid(divm,
                                                                  vv,
                                                                  out2);

      if (!out2.WriteToFile("1d_accuracy_out2.dat", "x", "~p'(x)")) {
        std::cerr << "Comp Grad Pol could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::UniStgGrid1D vv_prime(aa, bb, nn);

      vv_prime.BindScalarField(PolynomialDerivative1, order_in);

      if (!vv_prime.WriteToFile("1d_accuracy_vv_prime.dat", "x", "p'(x)")) {
        std::cerr << "Div Poly could not be written to file." << std::endl;
        return EXIT_FAILURE;
      }

      mtk::Real relative_norm_2_error{};

      relative_norm_2_error =
        mtk::BLASAdapter::RelNorm2Error(out2.discrete_field(),
                                        vv_prime.discrete_field(),
                                        vv_prime.num_cells_x());

      output_dat_file <<  " & " << relative_norm_2_error;
    }

    output_dat_file <<  "\\\\" << std::endl;
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
