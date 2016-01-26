# Makefile for example: divergence_operators_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

divergence_operators_1d.exe: \
  divergence_operators_1d.cc \
  divergence_operators_1d.o
	$(LINKER) divergence_operators_1d.o $(LIBS) \
  -o divergence_operators_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

divergence_operators_1d.o: divergence_operators_1d.cc $(MTK_LIB)
	$(CC) divergence_operators_1d.cc $(CCFLAGS) -o divergence_operators_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) divergence_operators_1d.exe >& \
          divergence_operators_1d.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
