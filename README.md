# The Mimetic Methods Toolkit (MTK)

By: **Eduardo J. Sanchez, PhD - esanchez at mail dot sdsu dot edu**

## 1. Description

We define numerical methods that are based on discretizations preserving the
properties of their continuous counterparts to be **mimetic**.

The **Mimetic Methods Toolkit (MTK)** is a C++11 library for mimetic numerical
methods. It is a set of classes for **mimetic interpolation**, **mimetic
quadratures**, and **mimetic finite difference** methods for the **numerical
solution of ordinary and partial differential equations**.

## 2. Dependencies

This README file assumes all of these dependencies are installed in the
following folder:

```
$(HOME)/Libraries/
```

In this version, the MTK optionally uses ATLAS-optimized BLAS and LAPACK
routines for the internal computation on some of the layers. However, ATLAS
requires both BLAS and LAPACK in order to create their optimized distributions.
Therefore, the following dependencies tree arises:

### For Linux:

1. LAPACK - Available from: http://www.netlib.org/lapack/
  1. BLAS - Available from: http://www.netlib.org/blas/

2. GLPK - Available from: https://www.gnu.org/software/glpk/

3. (Optional) ATLAS - Available from: http://math-atlas.sourceforge.net/
  1. LAPACK - Available from: http://www.netlib.org/lapack/
    1. BLAS - Available from: http://www.netlib.org/blas

4. (Optional) Valgrind - Available from: http://valgrind.org/

5. (Optional) Doxygen - Available from http://www.stack.nl/~dimitri/doxygen/

### For OSX:

1. GLPK - Available from: https://www.gnu.org/software/glpk/

## 3. Installation

### PART 1. CONFIGURATION OF THE MAKEFILE.

The following steps are required to build and test the MTK. Please use the
accompanying `Makefile.inc` file, which should provide a solid template to
start with. The following command provides help on the options for make:

```
$ make help
-----
Makefile for the MTK.

Options are:
- all: builds the library, the tests, and examples.
- mtklib: builds the library.
- test: builds the test files.
- example: builds the examples.

- testall: runs all the tests.

- gendoc: generates the documentation for the library.

- clean: cleans all the generated files.
- cleanlib: cleans the generated archive and object files.
- cleantest: cleans the generated tests executables.
- cleanexample: cleans the generated examples executables.
-----
```

### PART 2. BUILD THE LIBRARY.

```
$ make
```

If successful you'll read (before building the tests and examples):
```
----- Library created! Check in /home/ejspeiro/Dropbox/MTK/lib
```

## 4. Contact, Support, and Credits

The GitHub repository is: https://github.com/ejspeiro/MTK

The MTK is developed by researchers and adjuncts to the
[Computational Science Research Center (CSRC)](http://www.csrc.sdsu.edu/)
at [San Diego State University (SDSU)](http://www.sdsu.edu/).

Currently the developers are:

- **Eduardo J. Sanchez, PhD - esanchez at mail dot sdsu dot edu** - @ejspeiro
- Jose E. Castillo, PhD - jcastillo at mail dot sdsu dot edu
- Guillermo F. Miranda, PhD - unigrav at hotmail dot com

### 4.1. Acknowledgements and Contributions

The authors would like to acknowledge valuable advising, feedback,
and actual contributions from research personnel at the Computational Science
Research Center (CSRC) at San Diego State University (SDSU). Their input was
important to the fruition of this work. Specifically, our thanks go to
(alphabetical order):

- Mohammad Abouali, PhD
- Dany De Cecchis, PhD
- Otilio Rojas, PhD
- Julia Rossi, PhD
- Christopher P. Paolini, PhD - paolini at engineering dot sdsu dot edu
- Johnny Corbino.
- Raul Vargas-Navarro.

## 5. Referencing This Work

Please reference this work as follows:
```
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
```

Finally, please feel free to contact me with suggestions or corrections:

**Eduardo J. Sanchez, PhD - esanchez at mail dot sdsu dot edu** - @ejspeiro

Thanks and happy coding!
