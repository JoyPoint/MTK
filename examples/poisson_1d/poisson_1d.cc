/*!
\file poisson_1d.cc

\brief Poisson Equation on a 1D Uniform Staggered Grid with Robin BCs.

We solve:

\f[
\nabla^2p(x) = -s(x),
\f]

for \f$ x \in \Omega = [a,b] = [0,1] \f$.

The source term function is defined as

\f[
s(x) = \frac{\lambda^2\exp(\lambda x)}{\exp(\lambda) - 1}
\f]

where \f$ \lambda = -1 \f$ is a parameter.

We consider Robin's boundary conditions of the form:

\f[
\alpha p(a) - \beta p'(a) = \omega,
\f]

\f[
\alpha p(b) + \beta p'(b) = \epsilon.
\f]

The analytical solution for this problem is given by

\f[
p(x) = \frac{\exp(\lambda x) - 1}{\exp(\lambda) - 1}.
\f]

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\author: Raul Vargas--Navarro - vargasna at rohan dot sdsu dot edu
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

#if __cplusplus == 201103L

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>

#include "mtk.h"

mtk::Real Source(mtk::Real xx) {

  mtk::Real lambda = -1.0;

  return lambda*lambda*exp(lambda*xx)/(exp(lambda) - 1.0);
}

mtk::Real KnownSolution(mtk::Real xx) {

  mtk::Real lambda = -1.0;

  return (exp(lambda*xx) - 1.0)/(exp(lambda) - 1.0);
}

int main () {

  std::cout << "Example: Poisson Equation on a 1D Uniform Staggered Grid ";
  std::cout << "with Robin BCs." << std::endl;

  /// 1. Describe the problem of interest.

  mtk::Real lambda = -1.0;
  mtk::Real alpha = -exp(lambda);
  mtk::Real beta = (exp(lambda) - 1.0)/lambda;
  mtk::Real omega = -1.0;
  mtk::Real epsilon = 0.0;

  /// 2. Discretize space.

  mtk::Real west_bndy_x = 0.0;
  mtk::Real east_bndy_x = 1.0;
  int num_cells_x = 5;

  mtk::UniStgGrid1D comp_sol(west_bndy_x, east_bndy_x, num_cells_x);

  /// 3. Create mimetic operators.

  int order_of_accuracy{2};  // Desired order of accuracy for approximation.

  mtk::Grad1D grad;  // Mimetic gradient operator.

  mtk::Lap1D lap;  // Mimetic Laplacian operator.

  if (!lap.ConstructLap1D(order_of_accuracy)) {
    std::cerr << "Mimetic lap could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::DenseMatrix lapm(lap.ReturnAsDenseMatrix(comp_sol));

  std::cout << "Mimetic Laplacian operator: " << std::endl;
  std::cout << lapm << std::endl;

  if (!grad.ConstructGrad1D(order_of_accuracy)) {
    std::cerr << "Mimetic grad could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::DenseMatrix gradm(grad.ReturnAsDenseMatrix(comp_sol));

  std::cout << "Mimetic gradient operator: " << std::endl;
  std::cout << gradm << std::endl;

  /// 4. Create grid for source term.

  mtk::UniStgGrid1D source(west_bndy_x, east_bndy_x, num_cells_x);

  source.BindScalarField(Source);

  std::cout << source << std::endl;

  /// 5. Apply Boundary Conditions to both operator and source term.

  // Since we need to approximate the first derivative times beta, we must use
  // the approximation of the gradient at the boundary. We could extract them
  // from the gradient operator as packed in the grad object. BUT, since we have
  // generated at matrix containing this operator, we can extract these from the
  // matrix.

  // Array containing the coefficients for the west boundary condition.
  std::vector<mtk::Real> west_coeffs;

  for (auto ii = 0; ii < grad.num_bndy_coeffs(); ++ii) {
    west_coeffs.push_back(-beta*gradm.GetValue(0, ii));
  }

  // Array containing the coefficients for the east boundary condition.
  std::vector<mtk::Real> east_coeffs;

  for (auto ii = 0; ii < grad.num_bndy_coeffs(); ++ii) {
    east_coeffs.push_back(beta*gradm.GetValue(gradm.num_rows() - 1,
                                              gradm.num_cols() - 1 - ii));
  }

  // To impose the Dirichlet condition, we simple add its coefficient to the
  // first entry of the west, and the last entry of the east array.

  west_coeffs[0] += alpha;

  east_coeffs[0] += alpha;

  // Now that we have the coefficients that should be in the operator, we create
  // a boundary condition descriptor object, which will encapsulate the
  // complexity of assigning them in the matrix, to complete the construction of
  // the mimetic operator.

  mtk::BCDesc1D::ImposeOnOperator(lapm, west_coeffs, east_coeffs);

  std::cout << "Mimetic Laplacian with Robin conditions:" << std::endl;
  std::cout << lapm << std::endl;

  mtk::BCDesc1D::ImposeOnGrid(source, omega, epsilon);

  std::cout << "Source term with imposed BCs:" << std::endl;
  std::cout << source << std::endl;

  source.WriteToFile("poisson_1d_source.dat", "x", "s(x)");

  /// 6. Solve the problem.

  int info{mtk::LAPACKAdapter::SolveDenseSystem(lapm, source)};

  if (!info) {
    std::cout << "System solved! Problem solved!" << std::endl;
    std::cout << std::endl;
  }
  else {
    std::cerr << "Something wrong solving system! info = " << info << std::endl;
    std::cerr << "Exiting..." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Computed solution:" << std::endl;
  std::cout << source << std::endl;

  source.WriteToFile("poisson_1d_comp_sol.dat", "x", "~u(x)");

  /// 7. Compare computed solution against known solution.

  mtk::UniStgGrid1D known_sol(west_bndy_x, east_bndy_x, num_cells_x);

  known_sol.BindScalarField(KnownSolution);

  std::cout << "known_sol =" << std::endl;
  std::cout << known_sol << std::endl;

  known_sol.WriteToFile("poisson_1d_known_sol.dat", "x", "u(x)");

  mtk::Real relative_norm_2_error{};  // Relative norm 2 of the error.

  relative_norm_2_error =
    mtk::BLASAdapter::RelNorm2Error(source.discrete_field_u(),
                                    known_sol.discrete_field_u(),
                                    known_sol.num_cells_x());

  std::cout << "relative_norm_2_error = ";
  std::cout << relative_norm_2_error << std::endl;
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
