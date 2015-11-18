# Makefile for the tests folder.

include ../Makefile.inc

all: mtk_dense_matrix_test \
	mtk_blas_adapter_test mtk_lapack_adapter_test mtk_glpk_adapter_test \
	mtk_uni_stg_grid_1d_test mtk_uni_stg_grid_2d_test \
	mtk_grad_1d_test mtk_grad_2d_test \
	mtk_div_1d_test mtk_div_2d_test \
	mtk_lap_1d_test mtk_lap_2d_test \
	mtk_bc_descriptor_2d_test \
	mtk_interp_1d_test

EXE = \
	mtk_dense_matrix_test.exe \
	mtk_blas_adapter_test.exe \
	mtk_lapack_adapter_test.exe \
	mtk_glpk_adapter_test.exe \
	mtk_uni_stg_grid_1d_test.exe \
	mtk_uni_stg_grid_2d_test.exe \
	mtk_grad_1d_test.exe \
	mtk_grad_2d_test.exe \
	mtk_div_1d_test.exe \
	mtk_div_2d_test.exe \
	mtk_lap_1d_test.exe \
	mtk_lap_2d_test.exe \
	mtk_bc_descriptor_2d_test.exe \
	mtk_interp_1d_test.exe

#   Compiling statements.
#   __________________________________________________________________

%.o: %.cc $(MTK_LIB)
	$(CC) $(CCFLAGS) $<

#   Linking statements.
#   __________________________________________________________________

mtk_dense_matrix_test: mtk_dense_matrix_test.o mtk_dense_matrix_test.cc
	$(LINKER) mtk_dense_matrix_test.o $(LIBS) -o mtk_dense_matrix_test.exe
	@echo ' '

mtk_blas_adapter_test: mtk_blas_adapter_test.o mtk_blas_adapter_test.cc
	$(LINKER) mtk_blas_adapter_test.o $(LIBS) -o mtk_blas_adapter_test.exe
	@echo ' '

mtk_lapack_adapter_test: mtk_lapack_adapter_test.o mtk_lapack_adapter_test.cc
	$(LINKER) mtk_lapack_adapter_test.o $(LIBS) -o mtk_lapack_adapter_test.exe
	@echo ' '

mtk_glpk_adapter_test: mtk_glpk_adapter_test.o mtk_glpk_adapter_test.cc
	$(LINKER) mtk_glpk_adapter_test.o $(LIBS) -o mtk_glpk_adapter_test.exe
	@echo ' '

mtk_uni_stg_grid_1d_test: mtk_uni_stg_grid_1d_test.o mtk_uni_stg_grid_1d_test.cc
	$(LINKER) mtk_uni_stg_grid_1d_test.o $(LIBS) -o mtk_uni_stg_grid_1d_test.exe
	@echo ' '

mtk_uni_stg_grid_2d_test: mtk_uni_stg_grid_2d_test.o mtk_uni_stg_grid_2d_test.cc
	$(LINKER) mtk_uni_stg_grid_2d_test.o $(LIBS) -o mtk_uni_stg_grid_2d_test.exe
	@echo ' '

mtk_grad_1d_test: mtk_grad_1d_test.o mtk_grad_1d_test.cc
	$(LINKER) mtk_grad_1d_test.o $(LIBS) -o mtk_grad_1d_test.exe
	@echo ' '

mtk_grad_2d_test: mtk_grad_2d_test.o mtk_grad_2d_test.cc
	$(LINKER) mtk_grad_2d_test.o $(LIBS) -o mtk_grad_2d_test.exe
	@echo ' '

mtk_div_1d_test: mtk_div_1d_test.o mtk_div_1d_test.cc
	$(LINKER) mtk_div_1d_test.o $(LIBS) -o mtk_div_1d_test.exe
	@echo ' '

mtk_div_2d_test: mtk_div_2d_test.o mtk_div_2d_test.cc
	$(LINKER) mtk_div_2d_test.o $(LIBS) -o mtk_div_2d_test.exe
	@echo ' '

mtk_lap_1d_test: mtk_lap_1d_test.o mtk_lap_1d_test.cc
	$(LINKER) mtk_lap_1d_test.o $(LIBS) -o mtk_lap_1d_test.exe
	@echo ' '

mtk_lap_2d_test: mtk_lap_2d_test.o mtk_lap_2d_test.cc
	$(LINKER) mtk_lap_2d_test.o $(LIBS) -o mtk_lap_2d_test.exe
	@echo ' '

mtk_bc_descriptor_2d_test: mtk_bc_descriptor_2d_test.o \
  mtk_bc_descriptor_2d_test.cc
	$(LINKER) mtk_bc_descriptor_2d_test.o $(LIBS) \
	  -o mtk_bc_descriptor_2d_test.exe
	@echo ' '

mtk_interp_1d_test: mtk_interp_1d_test.o mtk_interp_1d_test.cc
	$(LINKER) mtk_interp_1d_test.o $(LIBS) -o mtk_interp_1d_test.exe
	@echo ' '

#   Utilities.
#   __________________________________________________________________

testall: mtk_dense_matrix_test_exe \
	\
	mtk_blas_adapter_test_exe \
	mtk_lapack_adapter_test_exe \
	mtk_glpk_adapter_test_exe \
	\
	mtk_uni_stg_grid_1d_test_exe \
	mtk_uni_stg_grid_2d_test_exe \
	\
	mtk_grad_1d_test_exe \
	mtk_grad_2d_test_exe \
	mtk_div_1d_test_exe \
	mtk_div_2d_test_exe \
	mtk_lap_1d_test_exe \
	mtk_lap_2d_test_exe \
	mtk_bc_descriptor_2d_test_exe \
	mtk_interp_1d_test_exe

mtk_dense_matrix_test_exe:
	./mtk_dense_matrix_test.exe

mtk_blas_adapter_test_exe:
	./mtk_blas_adapter_test.exe

mtk_lapack_adapter_test_exe:
	./mtk_lapack_adapter_test.exe

mtk_glpk_adapter_test_exe:
	./mtk_glpk_adapter_test.exe

mtk_uni_stg_grid_1d_test_exe:
	./mtk_uni_stg_grid_1d_test.exe

mtk_uni_stg_grid_2d_test_exe:
	./mtk_uni_stg_grid_2d_test.exe

mtk_grad_1d_test_exe:
	./mtk_grad_1d_test.exe

mtk_grad_2d_test_exe:
	./mtk_grad_2d_test.exe

mtk_div_1d_test_exe:
	./mtk_div_1d_test.exe

mtk_div_2d_test_exe:
	./mtk_div_2d_test.exe

mtk_lap_1d_test_exe:
	./mtk_lap_1d_test.exe

mtk_lap_2d_test_exe:
	./mtk_lap_2d_test.exe

mtk_bc_descriptor_2d_test_exe:
	./mtk_bc_descriptor_2d_test.exe

mtk_interp_1d_test_exe:
	./mtk_interp_1d_test.exe

memcheck:
	valgrind $(MEMCHECK_OPTS) \
	./mtk_dense_matrix_test.exe >& mtk_dense_matrix_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_blas_adapter_test.exe >& mtk_blas_adapter_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_lapack_adapter_test.exe >& mtk_lapack_adapter_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_glpk_adapter_test.exe >& mtk_glpk_adapter_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_uni_stg_grid_1d_test.exe >& mtk_uni_stg_grid_1d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_uni_stg_grid_2d_test.exe >& mtk_uni_stg_grid_2d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_grad_1d_test.exe >& mtk_grad_1d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_grad_2d_test.exe >& mtk_grad_2d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_div_1d_test.exe >& mtk_div_1d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_div_2d_test.exe >& mtk_div_2d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_lap_1d_test.exe >& mtk_lap_1d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_lap_2d_test.exe >& mtk_lap_2d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_bc_descriptor_2d_test.exe >& mtk_bc_descriptor_2d_test.memcheck
	valgrind $(MEMCHECK_OPTS) \
	./mtk_interp_1d_test.exe >& mtk_interp_1d_test.memcheck

clean:
	rm -f *~ *.o *.mps *.dat *.png *.memcheck $(EXE)
