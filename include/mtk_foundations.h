/*!
\file mtk_foundations.h

\brief Fundamental definitions to be used across all classes of the MTK.

This file contains the fundamental definitions that classes of the MTK rely on
to be implemented. Examples of these definitions are the definition of
fundamental data types, and global variables affecting the construction of
mimetic operators, among others.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\todo Test selective precision mechanisms.
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

#ifndef MTK_INCLUDE_ROOTS_H_
#define MTK_INCLUDE_ROOTS_H_

/*!
\namespace mtk

\brief Mimetic Methods Toolkit namespace.
*/
namespace mtk {

/*!
\def MTK_PRECISION_DOUBLE

\ingroup c01-roots

\brief Defined from Makefile, decides between a single- or double-precision MTK.

\warning Not defined anywhere in the source code, but on file Makefile.inc.
*/

/*!
\typedef Real

\ingroup c01-roots

\brief Users can simply change this to build a double- or single-precision MTK.

\warning Defined as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.
*/
#ifdef MTK_PRECISION_DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

/*!
\var kZero

\ingroup c01-roots

\brief MTK's zero defined according to selective compilation.

\warning Declared as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.

\var kOne

\brief MTK's one defined according to selective compilation.

\ingroup c01-roots

\warning Declared as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.

\var kTwo

\brief MTK's two defined according to selective compilation.

\ingroup c01-roots

\warning Declared as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.
*/
#ifdef MTK_PRECISION_DOUBLE
const double kZero{0.0};
const double kOne{1.0};
const double kTwo{2.0};
#else
const float kZero{0.0f};
const float kOne{1.0f};
const float kTwo{2.0f};
#endif

/*!
\var kDefaultTolerance

\ingroup c01-roots

\brief Considered tolerance for comparisons in numerical methods.

\warning Declared as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.
*/
#ifdef MTK_PRECISION_DOUBLE
const double kDefaultTolerance{1e-7};
#else
const float kDefaultTolerance{1e-7f};
#endif

/*!
\var kDefaultMimeticThreshold

\ingroup c01-roots

\brief Default tolerance for higher-order mimetic operators.

\warning Declared as double if MTK_PRECISION_DOUBLE is defined on Makefile.inc.
*/
#ifdef MTK_PRECISION_DOUBLE
const double kDefaultMimeticThreshold{1e-6};
#else
const float kDefaultMimeticThreshold{1e-6f};
#endif

/*!
\var kDefaultOrderAccuracy

\ingroup c01-roots

\brief Default order of accuracy for mimetic operators.
*/
const int kDefaultOrderAccuracy{2};

/*!
\var kCriticalOrderAccuracyGrad

\ingroup c01-roots

\brief At this order (and higher) we must use the CBSA to construct gradients.
*/
const int kCriticalOrderAccuracyGrad{10};

/*!
\var kCriticalOrderAccuracyDiv

\ingroup c01-roots

\brief At this order (and higher) we must use the CBSA to construct divergences.
*/
const int kCriticalOrderAccuracyDiv{8};
}
#endif  // End of: MTK_INCLUDE_ROOTS_H_
