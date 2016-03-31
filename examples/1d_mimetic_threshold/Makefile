# Makefile for example: 1d_mimetic_threshold.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_mimetic_threshold.exe: 1d_mimetic_threshold.cc 1d_mimetic_threshold.o
	$(LINKER) 1d_mimetic_threshold.o $(LIBS) -o 1d_mimetic_threshold.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_mimetic_threshold.o: 1d_mimetic_threshold.cc $(MTK_LIB)
	$(CC) 1d_mimetic_threshold.cc $(CCFLAGS) -o 1d_mimetic_threshold.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_mimetic_threshold.exe >& \
  1d_mimetic_threshold.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
