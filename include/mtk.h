/*!
\file mtk.h

\brief Includes the entire API.

This file contains every required header file, thus containing the entire API.
In this way, client codes only have to instruct #include "mtk.h".

\warning IT IS EXTREMELY IMPORTANT THAT THE HEADERS ARE ADDED TO THIS FILE IN A
SPECIFIC ORDER; THAT IS, CONSIDERING THE DEPENDENCE BETWEEN THE CLASSES THESE
CONTAIN!

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at mail dot sdsu dot edu
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
/*!
\mainpage Introduction

We define numerical methods that are based on discretizations preserving the
properties of their continuum counterparts to be **mimetic**.

The **Mimetic Methods Toolkit (MTK)** is a C++ library for mimetic numerical
methods. It is a set of classes for **mimetic quadratures**, **mimetic
interpolation**, and **mimetic discretization methods** for the numerical
solution of ordinary and partial differential equations.

\section section_mtk_concerns MTK Concerns

Since collaborative development efforts are definitely important in achieving
the level of generality we intend the library to possess, we have divided the
library's source code according to the designated purpose the classes possess
within the library. These divisions (or concerns) are grouped by layers, and are
hierarchically related by the dependence they have among them.

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

\section section_flavors MTK Flavors

The MTK collection of wrappers is:

-# MMTK: MATLAB wrappers collection for MTK; intended for sequential
computations.

Others are being designed and developed.

\section section_authors Contact, Support and Credits

The MTK is developed by researchers and adjuncts to the
[Computational Science Research Center (CSRC)](http://www.csrc.sdsu.edu/)
at [San Diego State University (SDSU)](http://www.sdsu.edu/).

Developers are members of:

1. Mimetic Numerical Methods Research and Development Group.
2. Computational Geoscience Research and Development Group.
3. Ocean Modeling Research and Development Group.

Currently the developers are:

-# **Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu** - ejspeiro
-# Jose E. Castillo, Ph.D. - jcastillo at mail dot sdsu dot edu
-# Guillermo F. Miranda, Ph.D. - unigrav at hotmail dot com
-# Christopher P. Paolini, Ph.D. - paolini at engineering dot sdsu dot edu

-# Angel Boada.
-# Johnny Corbino.
-# Raul Vargas-Navarro.

\section subsection_acknowledgements Acknowledgements and Contributions

The authors would like to acknowledge valuable advising, contributions and
feedback, from research personnel at the Computational Science Research Center
at San Diego State University, which were vital to the fruition of this work.
Specifically, our thanks go to (alphabetical order):

-# Mohammad Abouali, Ph.D.
-# Dany De Cecchis, Ph.D.
-# Julia Rossi.

\page section_prog_tools Programming Tools

The development of MTK has been made possible through the use of the following
applications:
-# Editor: Kate - KDE Advanced Text Editor. Version 3.13.3. Using KDE
Development Platform 4.13.3 (C) 2000-2005 The Kate Authors.
-# Compiler: gcc version 4.4.5 (Ubuntu/Linaro 4.4.4-14ubuntu5). Copyright (C)
2013 Free Software Foundation, Inc.
-# Debugger: GNU gdb (Ubuntu 7.7.1-0ubuntu5~14.04.2) 7.7.1. Copyright (C) 2014
Free Software Foundation, Inc.

\page section_license_mod Licensing and Modifications

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

\page page_readme Read Me File and Installation Instructions

<pre>
# README File for the Mimetic Methods Toolkit (MTK)

By: **Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu**
    __________________________________________________________________

## 1. Description

We define numerical methods that are based on discretizations preserving the
properties of their continuum counterparts to be **mimetic**.

The **Mimetic Methods Toolkit (MTK)** is a C++ library for mimetic numerical
methods. It is arranged as a set of classes for **mimetic quadratures**,
**mimetic interpolation**, and **mimetic discretization** methods for the
numerical solution of ordinary and partial differential equations.
    __________________________________________________________________

## 2. Dependencies

This README assumes all of these dependencies are installed in the following
folder:

\verbatim
$(HOME)/Libraries/
\endverbatim

In this version, the MTK optionally uses ATLAS-optimized BLAS and LAPACK
routines for the internal computation on some of the layers. However, ATLAS
requires both BLAS and LAPACK in order to create their optimized distributions.
Therefore, the following dependencies tree arises:

### For Linux:

1. LAPACK - Available from: http://www.netlib.org/lapack/
  1. BLAS - Available from: http://www.netlib.org/blas/

2. (Optional) ATLAS - Available from: http://math-atlas.sourceforge.net/
  1. BLAS - Available from: http://www.netlib.org/blas/
  2. LAPACK - Available from: http://www.netlib.org/lapack/

3. (Optional) Valgrind - Available from: http://valgrind.org/

4. (Optional) Doxygen - Available from http://www.stack.nl/~dimitri/doxygen/

### For OS X:

There are no dependences for OS X.
    __________________________________________________________________

## 3. Installation

### PART 1. CONFIGURATION OF THE MAKEFILE.

The following steps are required the build and test the MTK. Please use the
accompanying Makefile.inc file, which should provide a solid template to start
with. The following command provides help on the options for make:

\verbatim
$ make help
-----
Makefile for the MTK.

Options are:
- make: builds only the library and the examples.
- all: builds the library, the examples and the documentation.
- mtklib: builds the library, i.e. generates the archive files.
- test: generates the tests.
- example: generates the examples.
- gendoc: generates the documentation for the library.

- clean: cleans ALL the generated files.
- cleanlib: cleans the generated archive and object files.
- cleantest: cleans the generated tests executables.
- cleanexample: cleans the generated examples executables.
-----
\endverbatim

### PART 2. BUILD THE LIBRARY.

\verbatim
$ make
\endverbatim

If successful you'll read (before building the examples):

\verbatim
----- Library created! Check in /home/ejspeiro/Dropbox/MTK/lib
\endverbatim

Examples and tests will also be built.
    __________________________________________________________________

## 4. Frequently Asked Questions

Q: Why haven't you guys implemented GBS to build the library?
A: I'm on it as we speak! ;)

Q: Is there any main reference when it comes to the theory on Mimetic Methods?
A: Yes! Check: http://www.csrc.sdsu.edu/mimetic-book

Q: Do I need to generate the documentation myself?
A: You can if you want to... but if you DO NOT want to, just go to our website.
    __________________________________________________________________

## 5. Contact, Support, and Credits

The MTK is developed by researchers and adjuncts to the
[Computational Science Research Center (CSRC)](http://www.csrc.sdsu.edu/)
at [San Diego State University (SDSU)](http://www.sdsu.edu/).

Developers are members of:

1. Mimetic Numerical Methods Research and Development Group.
2. Computational Geoscience Research and Development Group.
3. Ocean Modeling Research and Development Group.

Currently the developers are:

-# **Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu** - ejspeiro
-# Jose E. Castillo, Ph.D. - jcastillo at mail dot sdsu dot edu
-# Guillermo F. Miranda, Ph.D. - unigrav at hotmail dot com
-# Christopher P. Paolini, Ph.D. - paolini at engineering dot sdsu dot edu
-# Angel Boada.
-# Johnny Corbino.
-# Raul Vargas-Navarro.

Finally, please feel free to contact me with suggestions or corrections:

**Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu** - ejspeiro

Thanks and happy coding!
</pre>

\page page_architectures Tests and Test Architectures

Tests are given in the <a href="files.html">files list</a> section. They are
provided in the /tests/ folder within the distributed software.

In this page we intend to make a summary of all of the architectures in where
the MTK has been tested. The MTK is intended to be as portable as possible
throughout architectures. The following architectures have provided flawless
installations of the API and correct execution of the examples:

<pre>
1. Linux 3.2.0-23-generic-pae #36-Ubuntu SMP i386 GNU/Linux
   Intel(R) Pentium(R) M processor 1.73GHz 2048 KB of cache and stepping of 8
   gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)
</pre>

Further architectures will be tested!

\page page_examples Examples

Examples are given in the <a href="files.html">files list</a> section. They are
provided in the /examples/ folder within the distributed software.

\page page_ref_the User Manual, References and Theory

The main source of references for this work can be found in:

http://www.csrc.sdsu.edu/mimetic-book/

However, a .PDF copy of this manual can be found
<a href="../latex/refman.pdf">here</a>.
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
// #include "mtk_blas_adapter.h"
// #include "mtk_lapack_adapter.h"
// #include "mtk_glpk_adapter.h"

/*!
\defgroup c06-grids Grids.

\brief Uniform rectangular staggered grids.

Uniform rectangular staggered grids.
*/
// #include "mtk_uni_stg_grid_1d.h"
// #include "mtk_uni_stg_grid_2d.h"

/*!
\defgroup c07-mim_ops Mimetic operators.

\brief Mimetic operators.

Mimetic operators.
*/
// #include "mtk_grad_1d.h"
// #include "mtk_div_1d.h"
// #include "mtk_lap_1d.h"

#endif // End of: MTK_INCLUDE_MTK_H_
