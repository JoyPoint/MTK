# Makefile for example: 1d_poisson.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

1d_poisson.exe: 1d_poisson.cc 1d_poisson.o
	$(LINKER) 1d_poisson.o $(LIBS) -o 1d_poisson.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_poisson.o: 1d_poisson.cc $(MTK_LIB)
	$(CC) 1d_poisson.cc $(CCFLAGS) -o 1d_poisson.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_poisson.exe >& 1d_poisson.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
