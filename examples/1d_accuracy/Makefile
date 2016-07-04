# Makefile for example: 1d_accuracy.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_accuracy.exe: 1d_accuracy.cc 1d_accuracy.o
	$(LINKER) 1d_accuracy.o $(LIBS) -o 1d_accuracy.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_accuracy.o: 1d_accuracy.cc $(MTK_LIB)
	$(CC) 1d_accuracy.cc $(CCFLAGS) -o 1d_accuracy.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_accuracy.exe >& \
  1d_accuracy.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
