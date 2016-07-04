# Makefile for example: 1d_poisson_sensitivity_threshold_sensitivity_threshold.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_poisson_sensitivity_threshold.exe: 1d_poisson_sensitivity_threshold.cc \
	1d_poisson_sensitivity_threshold.o
	$(LINKER) 1d_poisson_sensitivity_threshold.o $(LIBS) \
	-o 1d_poisson_sensitivity_threshold.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_poisson_sensitivity_threshold.o: 1d_poisson_sensitivity_threshold.cc \
	$(MTK_LIB)
	$(CC) 1d_poisson_sensitivity_threshold.cc $(CCFLAGS) \
	-o 1d_poisson_sensitivity_threshold.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_poisson_sensitivity_threshold.exe \
	>& 1d_poisson_sensitivity_threshold.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
