# Makefile for example: minimalistic_poisson_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

minimalistic_poisson_1d.exe: minimalistic_poisson_1d.cc minimalistic_poisson_1d.o
	$(LINKER) minimalistic_poisson_1d.o $(LIBS) -o minimalistic_poisson_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

minimalistic_poisson_1d.o: minimalistic_poisson_1d.cc $(MTK_LIB)
	$(CC) minimalistic_poisson_1d.cc $(CCFLAGS) -o minimalistic_poisson_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) minimalistic_poisson_1d.exe >& minimalistic_poisson_1d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
