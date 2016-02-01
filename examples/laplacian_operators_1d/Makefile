# Makefile for example: laplacian_operators_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

laplacian_operators_1d.exe: \
  laplacian_operators_1d.cc \
  laplacian_operators_1d.o
	$(LINKER) laplacian_operators_1d.o $(LIBS) \
  -o laplacian_operators_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

laplacian_operators_1d.o: laplacian_operators_1d.cc $(MTK_LIB)
	$(CC) laplacian_operators_1d.cc $(CCFLAGS) -o laplacian_operators_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) laplacian_operators_1d.exe >& \
          laplacian_operators_1d.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
