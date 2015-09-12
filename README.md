# The Mimetic Methods Toolkit (MTK)

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
accompanying `Makefile.inc` file, which should provide a solid template to
start with. The following command provides help on the options for make:

```
$ make help
-----
Makefile for the MTK.

Options are:
- make: builds only the library and the examples.
- all: builds the library, the examples and the documentation.
- mtklib: builds the library, i.e. generates the archive files.
- tests: generates the tests.
- examples: generates the examples.
- gendoc: generates the documentation for the library.
- checkheaders: checks syntax of the header files.

- clean: cleans ALL the generated files.
- cleanlib: cleans the generated archive and object files.
- cleantests: cleans the generated tests executables.
- cleanexamples: cleans the generated examples executables.
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

Examples and tests will also be built.
    __________________________________________________________________

## 4. Frequently Asked Questions

Q: Why haven't you guys implemented GBS to build the library?
A: I'm on it as we speak! ;)

Q: When will the other flavors be ready?
A: Soon! I'm working on getting help on developing those.

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

- **Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu** - @ejspeiro
- Jose E. Castillo, Ph.D. - jcastillo at mail dot sdsu dot edu
- Guillermo F. Miranda, Ph.D. - unigrav at hotmail dot com
- Christopher P. Paolini, Ph.D. - paolini at engineering dot sdsu dot edu
- Angel Boada.
- Johnny Corbino.
- Raul Vargas-Navarro.

Finally, please feel free to contact me with suggestions or corrections:

**Eduardo J. Sanchez, Ph.D. - esanchez at mail dot sdsu dot edu** - @ejspeiro

Thanks and happy coding!
