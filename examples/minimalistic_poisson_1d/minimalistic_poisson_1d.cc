/*!
\file minimalistic_poisson_1d.cc

\brief Poisson Equation on a 1D Uniform Staggered Grid with Robin BCs.

We solve:

\f[
\nabla^2 p(x) = -s(x),
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

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "mtk.h"

mtk::Real Alpha(const mtk::Real &tt) {
  mtk::Real lambda = -1.0;
  return -exp(lambda);
}

mtk::Real Beta(const mtk::Real &tt) {
  mtk::Real lambda = -1.0;
  return (exp(lambda) - 1.0)/lambda;
};

mtk::Real Omega(const mtk::Real &tt) {
  return -1.0;
};

mtk::Real Epsilon(const mtk::Real &tt) {
  return 0.0;
};

mtk::Real Source(const mtk::Real &xx) {
  mtk::Real lambda = -1.0;
  return lambda*lambda*exp(lambda*xx)/(exp(lambda) - 1.0);
}

mtk::Real KnownSolution(const mtk::Real &xx) {
  mtk::Real lambda = -1.0;
  return (exp(lambda*xx) - 1.0)/(exp(lambda) - 1.0);
}

int main () {

  mtk::Real west_bndy_x = 0.0;
  mtk::Real east_bndy_x = 1.0;
  mtk::Real relative_norm_2_error{};
  int num_cells_x = 5;
  mtk::Grad1D grad;
  mtk::Lap1D lap;
  std::vector<mtk::Real> west_coeffs;
  std::vector<mtk::Real> east_coeffs;
  mtk::UniStgGrid1D source(west_bndy_x, east_bndy_x, num_cells_x);
  mtk::UniStgGrid1D comp_sol(west_bndy_x, east_bndy_x, num_cells_x);
  mtk::UniStgGrid1D known_sol(west_bndy_x, east_bndy_x, num_cells_x);
  if (!lap.ConstructLap1D()) {
    std::cerr << "Mimetic lap could not be built." << std::endl;
    return EXIT_FAILURE;
  }
  mtk::DenseMatrix lapm(lap.ReturnAsDenseMatrix(comp_sol));
  source.BindScalarField(Source);
  mtk::RobinBCDescriptor1D bcs;
  bcs.PushBackWestCoeff(Alpha);
  bcs.PushBackWestCoeff(Beta);
  bcs.PushBackEastCoeff(Alpha);
  bcs.PushBackEastCoeff(Beta);
  bcs.set_west_condition(Omega);
  bcs.set_east_condition(Epsilon);
  if (!bcs.ImposeOnLaplacianMatrix(lap, lapm)) {
    std::cerr << "BCs  could not be bound to the matrix." << std::endl;
    return EXIT_FAILURE;
  }
  bcs.ImposeOnGrid(source);
  int info{mtk::LAPACKAdapter::SolveDenseSystem(lapm, source)};
  if (info != 0) {
    std::cerr << "Something wrong solving system! info = " << info << std::endl;
    return EXIT_FAILURE;
  }

  source.WriteToFile("minimalistic_poisson_1d_comp_sol.dat", "x", "~u(x)");
  known_sol.BindScalarField(KnownSolution);
  relative_norm_2_error =
    mtk::BLASAdapter::RelNorm2Error(source.discrete_field(),
                                    known_sol.discrete_field(),
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
