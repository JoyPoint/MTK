/*!
\file poisson_2d.cc

\brief Poisson Equation on a 2D Uniform Staggered Grid with Robin BCs.

We solve:

\f[
\nabla^2 u(\mathbf{x}) = s(\mathbf{x}),
\f]

for \f$ \mathbf{x} \in \Omega = [0,1]^2 \f$.

The source term function is defined as

\f[
s(x,y) = xye^{-frac{1}{2}(x^2 + y^2)}(x^2 + y^2 - 6).
\f]

Let \f$ \Omega = S \cup N \cup W \cup E\f$. We consider Dirichlet and Neumann
boundary conditions of the following form:
\f[
\forall\mathbf{x}\in W: u(\mathbf{x}) = 0.
\f]

\f[
\forall\mathbf{x}\in E:
\f]

\f[
\forall\mathbf{x}\in S: u(\mathbf{x}) = 0.
\f]

\f[
\forall\mathbf{x}\in N:
\f]

The analytical solution for this problem is given by

\f[
u(x,y) = xye^{-frac{1}{2}(x^2 + y^2)}.
\f]

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

#include <iostream>
#include <fstream>
#include <cmath>

#include <vector>

#include "mtk.h"

mtk::Real Source(const mtk::Real &xx, const mtk::Real &yy) {

  mtk::Real x_squared{xx*xx};
  mtk::Real y_squared{yy*yy};
  mtk::Real aux{-0.5*(x_squared + y_squared)};

  return xx*yy*exp(aux)*(x_squared + y_squared - 6.0);
}

int main () {

  std::cout << "Example: Poisson Equation on a 2D Uniform Staggered Grid ";
  std::cout << "with Dirichlet and Neumann BCs." << std::endl;

  /// 1. Discretize space.
  mtk::Real west_bndy_x{0.0};
  mtk::Real east_bndy_x{1.0};
  mtk::Real south_bndy_y{0.0};
  mtk::Real north_bndy_y{1.0};
  int num_cells_x{5};
  int num_cells_y{5};

  mtk::UniStgGrid2D comp_sol(west_bndy_x, east_bndy_x, num_cells_x,
                             south_bndy_y, north_bndy_y, num_cells_y);

  /// 2. Create mimetic operator as a matrix.
  mtk::Lap2D lap;

  if (!lap.ConstructLap2D(comp_sol)) {
    std::cerr << "Mimetic Laplacian could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::DenseMatrix lapm(lap.ReturnAsDenseMatrix());

  if (!lapm.WriteToFile("poisson_2d_lapm.dat")) {
    std::cerr << "Laplacian matrix could not be written to disk." << std::endl;
    return EXIT_FAILURE;
  }

  /// 3. Create grid for source term.
  mtk::UniStgGrid2D source(west_bndy_x, east_bndy_x, num_cells_x,
                           south_bndy_y, north_bndy_y, num_cells_y);

  source.BindScalarField(Source);

  if(!source.WriteToFile("poisson_2d_source.dat", "x", "y", "s(x,y)")) {
    std::cerr << "Source term could not be written to disk." << std::endl;
    return EXIT_FAILURE;
  }
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
