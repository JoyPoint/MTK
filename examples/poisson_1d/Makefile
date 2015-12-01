# Makefile for example: poisson_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

poisson_1d.exe: poisson_1d.cc poisson_1d.o
	$(LINKER) poisson_1d.o $(LIBS) -o poisson_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

poisson_1d.o: poisson_1d.cc $(MTK_LIB)
	$(CC) poisson_1d.cc $(CCFLAGS) -o poisson_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) poisson_1d.exe >& poisson_1d.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
