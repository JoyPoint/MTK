# Makefile for example: 1d_poisson_supercritical.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_poisson_supercritical.exe: 1d_poisson_supercritical.cc \
	1d_poisson_supercritical.o
	$(LINKER) 1d_poisson_supercritical.o $(LIBS) -o 1d_poisson_supercritical.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_poisson_supercritical.o: 1d_poisson_supercritical.cc $(MTK_LIB)
	$(CC) 1d_poisson_supercritical.cc $(CCFLAGS) -o 1d_poisson_supercritical.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_poisson_supercritical.exe >&
		1d_poisson_supercritical.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
