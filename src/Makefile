# Makefile for the src folder.

include ../Makefile.inc

all: exetools datastructs nummethods grids mimops archive

EXETOOLS = \
	mtk_tools.o

DATASTRUCTS  = \
	mtk_matrix.o \
	mtk_dense_matrix.o

NUMMETHODS = \
	mtk_blas_adapter.o \
	mtk_lapack_adapter.o \
	mtk_glpk_adapter.o

GRIDS = \
	mtk_uni_stg_grid_1d.o \
	mtk_uni_stg_grid_2d.o

MIMOPS = \
	mtk_grad_1d.o \
	mtk_div_1d.o \
	mtk_lap_1d.o \
	mtk_bc_desc_1d.o \
	mtk_interp_1d.o \
	mtk_grad_2d.o \
	mtk_div_2d.o \
	mtk_lap_2d.o

exetools: $(EXETOOLS)

datastructs: $(DATASTRUCTS)

nummethods: $(NUMMETHODS)

grids: $(GRIDS)

mimops: $(MIMOPS)

#   Compiling statements.
#   __________________________________________________________________

%.o: %.cc
	$(CC) $(CCFLAGS) $<

#   Archiving statements.
#   __________________________________________________________________

# Create the library directory with "-p" so no prob if it exists:

archive: $(EXETOOLS) $(DATASTRUCTS) $(NUMMETHODS) $(GRIDS) $(MIMOPS)
	@(mkdir -p $(LIB))
	$(ARCH) $(ARCHFLAGS) $(MTK_LIB) \
	  $(EXETOOLS) $(DATASTRUCTS) $(NUMMETHODS) $(GRIDS) $(MIMOPS)
	$(RANLIB) $(MTK_LIB)
	@echo ' '
	@echo '----- Library created! Check for ' $(MTK_LIB)
	@echo ' '

#   Utilities.
#   __________________________________________________________________

clean:
	rm -f *~ *.o
