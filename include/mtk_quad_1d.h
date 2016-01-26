/*!
\file mtk_quad_1d.h

\brief Includes the definition of the class Quad1D.

Definition of a class that implements a 1D quadrature solver based on the
mimetic discretization of the gradient operator.

\sa mtk::Grad1D

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu

\todo Implement this class.
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

#ifndef MTK_INCLUDE_QUAD_1D_H_
#define MTK_INCLUDE_QUAD_1D_H_

#include <iostream>
#include <iomanip>

#include <vector>

namespace mtk {

/*!
\class Quad1D

\ingroup c07-mim_ops

\brief Implements a 1D mimetic quadrature.

This class implements a 1D quadrature solver based on the mimetic discretization
of the gradient operator.
*/
class Quad1D {
 public:
  /// \brief Output stream operator for printing.
  friend std::ostream& operator <<(std::ostream& stream, Quad1D &in);

  /// \brief Default constructor.
  Quad1D();

  /*!
  \brief Copy constructor.

  \param [in] div Given quadrature.
  */
  Quad1D(const Quad1D &quad);

  /// \brief Destructor.
  ~Quad1D();

  /*!
  \brief Get the degree of interpolating polynomial per sub-interval of domain.

  \return Degree of the interpolating polynomial per sub-interval of the domain.
  */
  int degree_approximation() const;

  /*!
  \brief Return collection of weights.

  \return Collection of weights.
  */
  Real *weights() const;

  /*!
  \brief Mimetic integration routine.

  \param [in] Integrand Real-valued function to integrate.
  \param [in] grid Given integration domain.

  \return Result of the integration.
  */
  Real Integrate(Real (*Integrand)(Real xx), UniStgGrid1D grid) const;

 private:
  int degree_approximation_;  ///< Degree of the interpolating polynomial.

  std::vector<Real> weights_;  ///< Collection of weights.
};
}
#endif  // End of: MTK_INCLUDE_QUAD_1D_H_
