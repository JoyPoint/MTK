# Makefile for example: divergence_operators_1d_mimetic_test.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

divergence_operators_1d_mimetic_test.exe: \
  divergence_operators_1d_mimetic_test.cc \
  divergence_operators_1d_mimetic_test.o
	$(LINKER) divergence_operators_1d_mimetic_test.o $(LIBS) \
  -o divergence_operators_1d_mimetic_test.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

divergence_operators_1d_mimetic_test.o: \
  divergence_operators_1d_mimetic_test.cc $(MTK_LIB)
	$(CC) divergence_operators_1d_mimetic_test.cc $(CCFLAGS) \
  -o divergence_operators_1d_mimetic_test.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) divergence_operators_1d_mimetic_test.exe >& \
          divergence_operators_1d_mimetic_test.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
