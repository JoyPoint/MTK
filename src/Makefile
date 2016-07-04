# Makefile for the src folder.

include ../Makefile.conf

all: exetools datastructs nummethods grids mimops archive

EXETOOLS = \
	mtk_tools.o
	
DATASTRUCTS  = \
	mtk_matrix.o \
	mtk_dense_matrix.o \
	mtk_crs_matrix.o

NUMMETHODS = \
	mtk_blas_adapter.o \
	mtk_lapack_adapter.o \
	mtk_glpk_adapter.o

GRIDS = \
	mtk_uni_stg_grid_1d.o \
	mtk_uni_stg_grid_2d.o \
	mtk_uni_stg_grid_3d.o

MIMOPS = \
	mtk_grad_1d.o \
	mtk_div_1d.o \
	mtk_lap_1d.o \
	mtk_robin_bc_descriptor_1d.o \
	mtk_interp_1d.o \
	mtk_grad_2d.o \
	mtk_div_2d.o \
	mtk_div_3d.o \
	mtk_curl_2d.o \
	mtk_lap_2d.o \
	mtk_lap_3d.o \
	mtk_robin_bc_descriptor_2d.o \
	mtk_grad_3d.o \
	mtk_operator_applicator.o

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
