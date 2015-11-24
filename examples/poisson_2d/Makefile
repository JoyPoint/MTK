# Makefile for example: poisson_2d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

poisson_2d.exe: poisson_2d.cc poisson_2d.o
	$(LINKER) poisson_2d.o $(LIBS) -o poisson_2d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

poisson_2d.o: poisson_2d.cc $(MTK_LIB)
	$(CC) poisson_2d.cc $(CCFLAGS) -o poisson_2d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) poisson_2d.exe >& poisson_2d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
