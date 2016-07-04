# Makefile for example: 1d_poisson_minimal.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_poisson_minimal.exe: 1d_poisson_minimal.cc 1d_poisson_minimal.o
	$(LINKER) 1d_poisson_minimal.o $(LIBS) -o 1d_poisson_minimal.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_poisson_minimal.o: 1d_poisson_minimal.cc $(MTK_LIB)
	$(CC) 1d_poisson_minimal.cc $(CCFLAGS) -o 1d_poisson_minimal.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_poisson_minimal.exe >& \
  1d_poisson_minimal.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
