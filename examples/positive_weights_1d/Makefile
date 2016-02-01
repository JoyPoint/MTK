# Makefile for example: positive_weights_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

positive_weights_1d.exe: positive_weights_1d.cc positive_weights_1d.o
	$(LINKER) positive_weights_1d.o $(LIBS) -o positive_weights_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

positive_weights_1d.o: positive_weights_1d.cc $(MTK_LIB)
	$(CC) positive_weights_1d.cc $(CCFLAGS) -o positive_weights_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) positive_weights_1d.exe >& \
          positive_weights_1d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
