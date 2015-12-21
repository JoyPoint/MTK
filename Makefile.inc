# Makefile setup file for the MTK.

SHELL := /bin/bash

# Please set the following variables up:

#   1. Absolute path to base directory of the MTK.
#   __________________________________________________________________

BASE = /home/esanchez/Dropbox/MTK

#   2. The machine (platform) identifier and required machine precision.
#   __________________________________________________________________

# Options are:
# - LINUX: A LINUX box installation.
# - OSX: Uses OS X optimized solvers.

PLAT = LINUX

# Options are:
# - SINGLE: Use 4 B floating point numbers.
# - DOUBLE: Use 8 B floating point numbers.

PRECISION = DOUBLE

#   3. Optimized solvers and operations by means of ATLAS in Linux?
#   __________________________________________________________________

# If you have selected OSX in step 1, then you don't need to worry about this.

# Options are ON xor OFF:

ATL_OPT = OFF

#   4. Paths to dependencies (header files for compiling).
#   __________________________________________________________________

# GLPK include path (soon to go):

GLPK_INC = $(HOME)/Libraries/glpk-4.35/include

# Linux: If ATLAS optimization is ON, users should only provide the path to
# ATLAS:

ATLAS_INC = $(HOME)/Libraries/ATLAS_3.8.4-CORE/include

# OS X: Do nothing.

#   5. Paths to dependencies (archive files for (static) linking).
#   __________________________________________________________________

# GLPK linking path (soon to go):

GLPK_LIB = $(HOME)/Libraries/glpk-4.35/lib/lib64/libglpk.a

# If optimization is OFF, then provide the paths for:

BLAS_LIB = $(HOME)/Libraries/BLAS-3.5.0/libblas.a
LAPACK_LIB = $(HOME)/Libraries/lapack-3.5.0/liblapack.a

# WARNING: Vendor libraries should be used whenever they are available.

# However, if optimization is ON, please provide the path the ATLAS' archive:

ATLAS_LIB   = $(HOME)/Libraries/ATLAS_3.8.4-CORE/ATLAS_3.8.4-BUILD-Citadel/lib

#   6. Compiler and its flags.
#   __________________________________________________________________

CC = g++

# Selective Verbose Execution for Quick Debugging. Options are defined per
# concern, and per data hierarchy on each concern.

# 0: NO verbose at all.

# 1: Enable verbose down to the 7th concern: messages.
# 2: Enable verbose down to the 7th concern: messages + scalar results.
# 3: Enable verbose down to the 7th concern. 1.1. + array results.
# 4: Enable verbose down to the 7th concern. 1.2. + matrix results.

# 5: Enable verbose down to the 6th concern: messages.
# 6: Enable verbose down to the 6th concern: messages + scalar results.
# 7: Enable verbose down to the 6th concern. 2.1. + array results.
# 8: Enable verbose down to the 6th concern. 2.2. + matrix results.

# 9: Enable verbose down to the 5th concern: messages.
# 10: Enable verbose down to the 5th concern: messages + scalar results.
# 11: Enable verbose down to the 5th concern. 3.1. + array results.
# 12: Enable verbose down to the 5th concern. 3.2. + matrix results.

# 13: Enable verbose down to the 4th concern: messages.
# 14: Enable verbose down to the 4th concern: messages + scalar results.
# 15: Enable verbose down to the 4th concern. 4.1. + array results.
# 16: Enable verbose down to the 4th concern. 4.2. + matrix results.

VERBOSE_LEVEL = 16

# Enable preventions. In the MTK, methods first validate their required
# pre-conditions in run-time. Similarly, in many points throughout the MTK
# codebase, different sanity checks are performed, as well. If this symbol is
# defined to be 0, the MTK will # perform no validations to enhance execution
# performance. Options are:
# - YES.
# - NO.

PERFORM_PREVENTIONS = YES

# Enables creation of LaTeX table verbosing the computation of mimetic weights.

VERBOSE_WEIGHTS = YES

# Flags recommended for release code:

CCFLAGS = -Wall -Werror -O3

# Flags recommended for debugging code:

CCFLAGS = -Wall -Werror -g

#   7. Archiver, its flags, and ranlib:
#   __________________________________________________________________

ARCH      = ar
ARCHFLAGS = cr

# If your system does not have "ranlib" then set: "RANLIB = echo":

RANLIB = echo

# But, if possible:

RANLIB = ranlib

#   8. Valgrind's memcheck options (optional):
#   __________________________________________________________________

MEMCHECK_OPTS = -v --tool=memcheck --leak-check=full --show-leak-kinds=all \
  --track-origins=yes --freelist-vol=20000000

# Done! User, please, do not mess with the definitions from this point on.

#   __________________________________________________________________
#   __________________________________________________________________
#   __________________________________________________________________

#   MTK-related.
#   __________________________________________________________________

SRC       = $(BASE)/src
INCLUDE   = $(BASE)/include
LIB       = $(BASE)/lib
MTK_LIB   = $(LIB)/libmtk.a
TESTS     = $(BASE)/tests
EXAMPLES  = $(BASE)/examples

#   Compiling-related.
#   __________________________________________________________________

CCFLAGS += -std=c++11 -fPIC \
  -DMTK_VERBOSE_LEVEL=$(VERBOSE_LEVEL) -I$(INCLUDE) -c

ifeq ($(PRECISION),DOUBLE)
  CCFLAGS += -DMTK_PRECISION_DOUBLE
else
  CCFLAGS += -DMTK_PRECISION_SINGLE
endif

ifeq ($(PERFORM_PREVENTIONS),YES)
  CCFLAGS += -DMTK_PERFORM_PREVENTIONS
endif

ifeq ($(VERBOSE_WEIGHTS),YES)
  CCFLAGS += -DMTK_VERBOSE_WEIGHTS
endif

# Only the GLPK is included because the other dependencies are coded in Fortran.

ifeq ($(ATL_OPT),ON)
  CCFLAGS  += -I$(GLPK_INC) $(ATLAS_INC)
else
  CCFLAGS  += -I$(GLPK_INC)
endif

#   Linking-related.
#   __________________________________________________________________

NOOPT_LIBS  = $(LAPACK_LIB) $(BLAS_LIB) -lm $(GLPK_LIB) -lstdc++

OPT_LIBS    = -L$(ATLAS_LIB) -latlas -llapack -lblas -lm -latlas -lstdc++

ifeq ($(PLAT),OSX)
  LINKER  = g++
  LINKER  += -framework Accelerate $(GLPK_LIB) $(MTK_LIB)
else
  ifeq ($(ATL_OPT),ON)
    LINKER  = g++
    LIBS = $(MTK_LIB)
    LIBS += $(OPT_LIBS)
  else
    LINKER  = gfortran
    LIBS = $(MTK_LIB)
    LIBS += $(NOOPT_LIBS)
  endif
endif

#   Documentation-related.
#   __________________________________________________________________

DOCGEN      = doxygen
DOCFILENAME = doc_config.dxcf
DOC         = $(BASE)/doc
DOCFILE     = $(BASE)/$(DOCFILENAME)
