# Makefile for example: 1d_positive_weights.exe

include ../../Makefile.conf

#   Linking statements.
#   __________________________________________________________________

1d_positive_weights.exe: 1d_positive_weights.cc 1d_positive_weights.o
	$(LINKER) 1d_positive_weights.o $(LIBS) -o 1d_positive_weights.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

1d_positive_weights.o: 1d_positive_weights.cc $(MTK_LIB)
	$(CC) 1d_positive_weights.cc $(CCFLAGS) -o 1d_positive_weights.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 1d_positive_weights.exe >& \
          1d_positive_weights.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe *.tex
