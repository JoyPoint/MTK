# Makefile setup file for MTK.

SHELL := /bin/bash

# Please set the following variables up:

#   1. Absolute path to base directory.
#   __________________________________________________________________

BASE = $(HOME)/Dropbox/MTK

#   2. The machine (platform) identifier and required precision.
#   __________________________________________________________________

# Options are:
# - LINUX: A LINUX box installation.
# - OSX: Soon!

PLAT = LINUX

# Options are:
# - SINGLE: Use 4 B floating point numbers.
# - DOUBLE: Use 8 B floating point numbers.

PRECISION = DOUBLE

#   3. Optimized solvers and operations by means of ATLAS in Linux?
#   __________________________________________________________________

# Options are ON xor OFF:

ATL_OPT = OFF

#   4. Paths to dependencies (header files for compiling).
#   __________________________________________________________________

# GLPK include path (soon to go):

GLPK_INC = $(HOME)/Libraries/glpk-4.55/include

# If ATLAS optimization is ON, users should only provide the path to ATLAS:

ATLAS_INC = $(HOME)/Libraries/ATLAS_3.8.4-CORE/include

#   5. Paths to dependencies (archive files for (static) linking).
#   __________________________________________________________________

# GLPK linking path (soon to go):

GLPK_LIB = $(HOME)/Libraries/glpk-4.55/lib/libglpk.a

# If optimization is OFF, then provide the paths for:

BLAS_LIB = $(HOME)/Libraries/BLAS/libblas.a
LAPACK_LIB = $(HOME)/Libraries/lapack-3.4.1/liblapack.a

# WARNING: Vendor libraries should be used whenever they are available.

# However, if optimization is ON, please provide the path the ATLAS' archive:

ATLAS_LIB   = $(HOME)/Libraries/ATLAS_3.8.4-CORE/ATLAS_3.8.4-BUILD-Citadel/lib

#   6. Compiler and its flags.
#   __________________________________________________________________

CC = colorgcc

# Debug Level. Options are:
# 0. NO debug at all NOR any run-time checks... be cautious!
# 1. Verbose (execution messages) AND run-time checks.
# 2. Level 1 plus intermediate scalar-valued results.
# 3. Level 2 plus intermediate array-valued results.

DEBUG_LEVEL = 3

# Flags recommended for release code:

CCFLAGS = -Wall -O2

# Flags recommended for debugging code:

CCFLAGS = -Wall -g

#   7. Archiver, its flags, and ranlib:
#   __________________________________________________________________

ARCH      = ar
ARCHFLAGS = cr

# If your system does not have "ranlib" then set: "RANLIB = echo":

RANLIB = echo

# But, if possible:

RANLIB = ranlib

#   8. Valgrind's memcheck options:
#   __________________________________________________________________

MEMCHECK_OPTS = -v --tool=memcheck --leak-check=full --show-leak-kinds=all \
--track-origins=yes --freelist-vol=20000000

# Done!

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

CCFLAGS += -std=c++11 -fPIC -DMTK_DEBUG_LEVEL=$(DEBUG_LEVEL) -I$(INCLUDE)  -c

ifeq ($(PRECISION),DOUBLE)
  CCFLAGS += -DMTK_PRECISION_DOUBLE
else
  CCFLAGS += -DMTK_PRECISION_SINGLE
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

ifeq ($(PLAT),MACOS)
  LINKER  = g++
  LINKER  += -framework Accelerate $(GLPK_LIB) $(MTK_LIB)
else
  ifeq ($(ATL_OPT),ON)
    LINKER  = g++
    LIBS = $(MTK_LIB)
    LIBS += $(OPT_LIBS)
  else
    LINKER  = g++
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
