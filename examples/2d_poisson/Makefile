# Makefile for example: 2d_poisson.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

2d_poisson.exe: 2d_poisson.cc 2d_poisson.o
	$(LINKER) 2d_poisson.o $(LIBS) -o 2d_poisson.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

2d_poisson.o: 2d_poisson.cc $(MTK_LIB)
	$(CC) 2d_poisson.cc $(CCFLAGS) -o 2d_poisson.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 2d_poisson.exe >& 2d_poisson.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
