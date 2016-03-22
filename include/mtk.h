/*!
\file mtk.h

\brief Includes the entire API.

This file contains every required header file, thus containing the entire API.
In this way, client codes only have to instruct #include "mtk.h".

\warning It is extremely important that the headers are added to this file in a
specific order; that is, considering the dependence between the classes these
contain.

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
/*!
\mainpage Introduction

We define numerical methods that are based on discretizations preserving the
properties of their continuous counterparts to be **mimetic**.

The **Mimetic Methods Toolkit (MTK)** is a C++11 library for mimetic numerical
methods. It is a set of classes for **mimetic interpolation**, **mimetic
quadratures**, and **mimetic finite difference** methods for the **numerical
solution of ordinary and partial differential equations**.

\section section_mtk_concerns MTK Concerns

Since collaborative development efforts are definitely important in achieving
the level of generality we intend the library to possess, we have divided the
library's source code according to the designated purpose the classes possess
within the library. These divisions (or **concerns**) are grouped by layers,
and are hierarchically related by the dependence they have among them.

One concern is said to depend on another one, if the classes the first concern
includes, rely on the classes the second concern includes.

In order of dependence these are:
-# Roots.
-# Enumerations.
-# Tools.
-# Data Structures.
-# Numerical Methods.
-# Grids.
-# Mimetic Operators.

\section section_wrappers MTK Wrappers

The MTK collection of wrappers is:

-# MMTK: MATLAB wrappers collection for MTK; intended for sequential
computations.

Others are being strongly considered.

\section section_authors Contact, Support and Credits

The GitHub repository is: https://github.com/ejspeiro/MTK

The MTK is developed by researchers and adjuncts to the
[Computational Science Research Center (CSRC)](http://www.csrc.sdsu.edu/)
at [San Diego State University (SDSU)](http://www.sdsu.edu/).

Currently the developers are:

- **Eduardo J. Sanchez, PhD - esanchez at mail dot sdsu dot edu** - ejspeiro
- Jose E. Castillo, PhD - jcastillo at mail dot sdsu dot edu
- Guillermo F. Miranda, PhD - unigrav at hotmail dot com
- Christopher P. Paolini, PhD - paolini at engineering dot sdsu dot edu
- Angel Boada.
- Johnny Corbino.
- Raul Vargas-Navarro.

\subsection subsection_acknowledgements Acknowledgements and Contributions

The authors would like to acknowledge valuable advising, contributions and
feedback, from research personnel at the Computational Science Research Center
at San Diego State University, which were vital to the fruition of this work.
Specifically, our thanks go to (alphabetical order):

-# Mohammad Abouali, Ph.D.
-# Dany De Cecchis, Ph.D.
-# Otilio Rojas, Ph.D.
-# Julia Rossi.

\page page_referencing Referencing This Work

Please reference this work as follows:
\verbatim
@article{Sanchez2014308,
  title = "The Mimetic Methods Toolkit: An object-oriented \{API\} for Mimetic
Finite Differences ",
  journal = "Journal of Computational and Applied Mathematics ",
  volume = "270",
  number = "",
  pages = "308 - 322",
  year = "2014",
  note = "Fourth International Conference on Finite Element Methods in
Engineering and Sciences (FEMTEC 2013) ",
  issn = "0377-0427",
  doi = "http://dx.doi.org/10.1016/j.cam.2013.12.046",
  url = "http://www.sciencedirect.com/science/article/pii/S037704271300719X",
  author = "Eduardo J. Sanchez and Christopher P. Paolini and Jose E. Castillo",
  keywords = "Object-oriented development",
  keywords = "Partial differential equations",
  keywords = "Application programming interfaces",
  keywords = "Mimetic Finite Differences "
}

@Inbook{Sanchez2015,
  author="Sanchez, Eduardo and Paolini, Christopher and Blomgren, Peter
and Castillo, Jose",
  editor="Kirby, M. Robert and Berzins, Martin and Hesthaven, S. Jan",
  chapter="Algorithms for Higher-Order Mimetic Operators",
  title="Spectral and High Order Methods for Partial Differential Equations
ICOSAHOM 2014: Selected papers from the ICOSAHOM conference, June 23-27, 2014,
Salt Lake City, Utah, USA",
  year="2015",
  publisher="Springer International Publishing",
  address="Cham",
  pages="425--434",
  isbn="978-3-319-19800-2",
  doi="10.1007/978-3-319-19800-2_39",
  url="http://dx.doi.org/10.1007/978-3-319-19800-2_39"
}
\endverbatim

\page page_readme Read Me File and Installation Instructions

\verbinclude ../README.md

\page section_prog_tools Programming Tools

The development of MTK has been made possible through the use of the following
applications:
-# Editor: Kate - KDE Advanced Text Editor. Version 3.13.3. Using KDE
Development Platform 4.13.3 (C) 2000-2005 The Kate Authors.
-# Debugger: GNU gdb (Ubuntu 7.7.1-0ubuntu5~14.04.2) 7.7.1. Copyright (C) 2014
Free Software Foundation, Inc.
-# Memory Profiler: valgrind-3.10.0.SVN.

See the section on test architectures for information about operating systems
and compilers used.

\page page_architectures Tests and Test Architectures

Tests are given in the <a href="files.html">files list</a> section. They are
provided in the /tests/ folder within the distributed software.

In this page we intend to make a summary of all of the architectures in where
the MTK has been tested. The MTK is intended to be as portable as possible
throughout architectures. The following architectures have provided flawless
installations of the API and correct execution of the tests and the examples:

\verbatim
1. Intel(R) Pentium(R) M CPU 1.73 GHz 2048 KB of cache and stepping of 8.
   Linux 3.2.0-23-generic-pae #36-Ubuntu SMP i386 GNU/Linux
   gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)

2. Intel(R) Core(TM) i7-4700MQ CPU 2.40 GHz 6144 KB of cache and stepping of 3.
   Linux 3.13.0-67-generic #110-Ubuntu SMP x86_64 GNU/Linux
   gcc version 4.8.4 (Ubuntu 4.4.4-2ubuntu1~14.04)

3. Intel(R) Core(TM) i7-4600U CPU 2.10 GHz 4096 KB of cache and a stepping of 1.
   Linux 3.16.7-29-desktop #1 SMP PREEMPT (6be6a97) x86_64 GNU/Linux
   openSUSE 13.2 (Harlequin) (x86_64)
   gcc (SUSE Linux) 4.8.3 20140627 [gcc-4_8-branch revision 212064]
\endverbatim

Further architectures will be tested!

\page page_ref_the User Manual, References and Theory

The main source of references for this work can be found in:

http://www.csrc.sdsu.edu/mimetic-book/

However, a .PDF copy of this manual can be found
<a href="../latex/refman.pdf">here</a>.

\page page_examples Examples

Examples are given in the <a href="files.html">files list</a> section. They are
provided in the /examples/ folder within the distributed software.

\page section_license_mod Licensing and Modifications

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

#ifndef MTK_INCLUDE_MTK_H_
#define MTK_INCLUDE_MTK_H_

/*!
\defgroup c01-roots Roots.

\brief Fundamental execution parameters and defined types.

Fundamental execution parameters and defined types.
*/
#include "mtk_roots.h"

/*!
\defgroup c02-enums Enumerations.

\brief Enumerations.

Enumerations.
*/
#include "mtk_enums.h"

/*!
\defgroup c03-execution_tools Execution tools.

\brief Tools to ensure execution correctness.

Tools to ensure execution correctness.
*/
#include "mtk_tools.h"

/*!
\defgroup c04-data_structures Data structures.

\brief Fundamental data structures.

Fundamental data structures.
*/
#include "mtk_matrix.h"
#include "mtk_dense_matrix.h"

/*!
\defgroup c05-num_methods Numerical methods.

\brief Adapter classes and auxiliary numerical methods.

Adapter classes and auxiliary numerical methods.
*/
#include "mtk_blas_adapter.h"
#include "mtk_lapack_adapter.h"
#include "mtk_glpk_adapter.h"

/*!
\defgroup c06-grids Grids.

\brief Uniform rectangular staggered grids.

Uniform rectangular staggered grids.
*/
#include "mtk_uni_stg_grid_1d.h"
#include "mtk_uni_stg_grid_2d.h"
#include "mtk_uni_stg_grid_3d.h"

/*!
\defgroup c07-mim_ops Mimetic operators.

\brief Mimetic operators.

Mimetic operators.
*/
#include "mtk_grad_1d.h"
#include "mtk_div_1d.h"
#include "mtk_lap_1d.h"
#include "mtk_robin_bc_descriptor_1d.h"
#include "mtk_quad_1d.h"
#include "mtk_interp_1d.h"

#include "mtk_grad_2d.h"
#include "mtk_div_2d.h"
#include "mtk_curl_2d.h"
#include "mtk_lap_2d.h"
#include "mtk_robin_bc_descriptor_2d.h"

#include "mtk_grad_3d.h"
#include "mtk_div_3d.h"
#include "mtk_lap_3d.h"
#include "mtk_robin_bc_descriptor_3d.h"

#include "mtk_operator_applicator.h"

#endif // End of: MTK_INCLUDE_MTK_H_
