/*!
\file mtk_lap_2d.cc

\brief Includes the implementation of the class Lap2D.

This class implements a 2D Laplacian operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm (CBSA).

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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_lap_2d.h"

mtk::Lap2D::Lap2D(): order_accuracy_(), mimetic_threshold_() {}

mtk::Lap2D::Lap2D(const Lap2D &lap):
  order_accuracy_(lap.order_accuracy_),
  mimetic_threshold_(lap.mimetic_threshold_) {}

mtk::Lap2D::~Lap2D() {}

bool mtk::Lap2D::ConstructLap2D(const mtk::UniStgGrid2D &grid,
                                int order_accuracy,
                                mtk::Real mimetic_threshold) {

  int NumCellsX{grid.num_cells_x()};
  int NumCellsY{grid.num_cells_y()};
  int aux{(NumCellsX + 2)*(NumCellsY + 2)};

  mtk::DenseMatrix gg;
  mtk::DenseMatrix dd;
  mtk::DenseMatrix lap;

   mtk::mtk::2DDiv *div = new mtk::mtk::2DDiv(kk_,mimetic_tol_);
   mtk::mtk::2DGrad *grad = new mtk::mtk::2DGrad(kk_,mimetic_tol_);

  


        East_ = East;
       West_=West;
      South_=South;
      North_=North;


   div = div->Construct2DDiv();
     if (div == nullptr) {
    return nullptr;
  }

  grad = grad->Construct2DGrad();
  if (grad == nullptr) {
    return nullptr;
  }


   gg=grad->ReturnAsMatrix(NumCellsX, West, East,NumCellsY, South, North);
   dd=div->ReturnAsMatrix(NumCellsX, West, East,NumCellsY, South, North);

   lap=new mtk::DenseMatrix(m,m);


  if (lap==nullptr)
   {
      cout << "Problem allocating memory for mtk::DenseMatrix"<<endl;
     return lap;
   }



   lap=lap->MatrixMultiplication(*dd,*gg);

  return lap;


}

mtk::DenseMatrix mtk::Lap2D::ReturnAsDenseMatrix() {

  return laplacian_;
}

mtk::Real* mtk::Lap2D::data() {

  return laplacian_.data();
}
