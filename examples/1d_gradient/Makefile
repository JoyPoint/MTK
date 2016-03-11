# Makefile for example: 1d_gradient.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

1d_gradient.exe: 1d_gradient.cc 1d_gradient.o
	$(LINKER) 1d_gradient.o $(LIBS) -o 1d_gradient.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_gradient.o: 1d_gradient.cc $(MTK_LIB)
	$(CC) 1d_gradient.cc $(CCFLAGS) -o 1d_gradient.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_gradient.exe >& 1d_gradient.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
