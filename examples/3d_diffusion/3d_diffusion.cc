/*!
\file 3d_diffusion.cc

\brief Diffusion Equation on a 3D Uniform Staggered Grid with Dirichlet BCs.

We solve:
\f[
\frac{\partial u}{\partial t} = \nabla^2 u(\mathbf{x}),
\f]
for \f$ \mathbf{x} \in \Omega = [0,1]^3 \f$.

We consider autonomous homogeneous Dirichlet boundary conditions.

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

mtk::Real InitialCondition(const mtk::Real &xx,
                           const mtk::Real &yy,
                           const mtk::Real &zz) {

  mtk::Real rr{0.3};

  mtk::Real aux{xx*xx + yy*yy + zz*zz};

  return (aux < rr? rr - aux: mtk::kZero);
}

int main () {

  std::cout << "Example: Diffusion Equation in 3D "
    "with Dirichlet BCs." << std::endl;

  /// 1. Discretize space.

  mtk::Real west_bndy_x{0.0};
  mtk::Real east_bndy_x{1.0};
  mtk::Real south_bndy_y{0.0};
  mtk::Real north_bndy_y{1.0};
  mtk::Real bottom_bndy_z{0.0};
  mtk::Real top_bndy_z{1.0};

  int num_cells_x{10};
  int num_cells_y{10};
  int num_cells_z{10};

  mtk::UniStgGrid3D comp_sol(west_bndy_x, east_bndy_x, num_cells_x,
                             south_bndy_y, north_bndy_y, num_cells_y,
                             bottom_bndy_z, top_bndy_z, num_cells_z);

  /// 2. Bind initial condition to the grid.

  comp_sol.BindScalarField(InitialCondition);

  if(!comp_sol.WriteToFile("3d_diffusion_comp_sol.dat",
                     "x",
                     "y",
                     "z",
                     "Initial u(x,y,z)")) {
    std::cerr << "Error writing to file." << std::endl;
    return EXIT_FAILURE;
  }

  /// 3. Create mimetic operator as a matrix.

  mtk::Lap3D lap;

  if (!lap.ConstructLap3D(comp_sol)) {
    std::cerr << "Mimetic Laplacian could not be built." << std::endl;
    return EXIT_FAILURE;
  }

  mtk::DenseMatrix lapm(lap.ReturnAsDenseMatrix());

  if (!lapm.WriteToFile("3d_diffusion_lapm.dat")) {
    std::cerr << "Laplacian matrix could not be written to disk." << std::endl;
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
