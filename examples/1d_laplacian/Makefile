# Makefile for example: 1d_laplacian.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_laplacian.exe: 1d_laplacian.cc 1d_laplacian.o
	$(LINKER) 1d_laplacian.o $(LIBS) -o 1d_laplacian.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_laplacian.o: 1d_laplacian.cc $(MTK_LIB)
	$(CC) 1d_laplacian.cc $(CCFLAGS) -o 1d_laplacian.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_laplacian.exe >& 1d_laplacian.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
