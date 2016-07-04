# Makefile for example: 1d_mimetic_measure.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_mimetic_measure.exe: 1d_mimetic_measure.cc 1d_mimetic_measure.o
	$(LINKER) 1d_mimetic_measure.o $(LIBS) -o 1d_mimetic_measure.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_mimetic_measure.o: 1d_mimetic_measure.cc $(MTK_LIB)
	$(CC) 1d_mimetic_measure.cc $(CCFLAGS) -o 1d_mimetic_measure.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_mimetic_measure.exe >& \
  1d_mimetic_measure.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
