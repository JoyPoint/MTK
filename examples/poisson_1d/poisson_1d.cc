/*!
\file poisson_1d.cc

\brief Poisson Equation on a 1D Uniform Staggered Grid with Robin BCs.

We solve:

\f[
-\nabla^2u(x) = s(x),
\f]

for \f$ x \in \Omega = [0,1] \f$.

The source term function is defined as

\f[
s(x) = \frac{-\lambda^2\exp(\lambda x)}{\exp(\lambda) - 1}
\f]

where \f$ \lambda = -1 \f$ is a parameter.

We consider Robin's boundary conditions of the form:

\f[
\alpha u(a) - \beta u'(a) = \textrm{west\_bndy\_cond\_val},
\f]

\f[
\alpha u(b) + \beta u'(b) = \textrm{east\_bndy\_cond\_val}.
\f]

The analytical solution for this problem is given by

\f[
u(x) = \frac{\exp(\lambda x) - 1}{\exp(\lambda) - 1}.
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

#include "mtk.h"

mtk::Real Source(mtk::Real xx) {

  mtk::Real lambda = -1.0;

  return -lambda*lambda*exp(lambda*xx)/(exp(lambda) - 1.0);
}

mtk::Real KnownSolution(mtk::Real xx) {

  mtk::Real lambda = -1.0;

  return (exp(lambda*xx) - 1.0)/(exp(lambda) - 1.0);
}

int main () {

  std::cout << "Example: Poisson Equation on a 1D Uniform Staggered Grid ";
  std::cout << "with Robin BCs." << std::endl;

  /// 1. Problem description.

  mtk::Real lambda = -1.0;
  mtk::Real alpha = -exp(lambda);
  mtk::Real beta = (exp(lambda) - 1.0)/lambda;
  mtk::Real west_bndy_value = -1.0;
  mtk::Real east_bndy_value = 0.0;

  /// 2. Spatial discretization.

  mtk::Real west_bndy_x = 0.0;
  mtk::Real east_bndy_x = 1.0;
  int num_cells_x = 5;

  mtk::UniStgGrid1D comp_sol(west_bndy_x, east_bndy_x, num_cells_x);

  /// 3. Mimetic Operators.

  int order_of_accuracy = 2;

  mtk::Grad1D grad; // Mimetic gradient operator.

  mtk::Lap1D lap; // Mimetic gradient operator.

  mtk::Real mimetic_threshold = 1e-6;

  if (!lap.ConstructLap1D(order_of_accuracy, mimetic_threshold)) {
    std::cerr << "Mimetic lap could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::DenseMatrix lapm(lap.ReturnAsDenseMatrix(comp_sol));

  std::cout << "Mimetic Laplacian operator: " << std::endl;
  std::cout << lapm << std::endl;

  if (!grad.ConstructGrad1D(order_of_accuracy, mimetic_threshold)) {
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

  source.WriteToFile("poisson_1d_source.dat", "x", "s(x)");

  /// 5. Apply Boundary Conditions to both operator and source term.

  /// 6. Solve the problem.

  /// 7. Compare computed solution against known solution.

  mtk::UniStgGrid1D known_sol(west_bndy_x, east_bndy_x, num_cells_x);

  known_sol.BindScalarField(KnownSolution);

  std::cout << known_sol << std::endl;

  known_sol.WriteToFile("poisson_1d_known_sol.dat", "x", "u(x)");
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
