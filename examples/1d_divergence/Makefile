# Makefile for example: 1d_divergence.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

1d_divergence.exe: 1d_divergence.cc 1d_divergence.o
	$(LINKER) 1d_divergence.o $(LIBS) -o 1d_divergence.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_divergence.o: 1d_divergence.cc $(MTK_LIB)
	$(CC) 1d_divergence.cc $(CCFLAGS) -o 1d_divergence.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_divergence.exe >& 1d_divergence.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
