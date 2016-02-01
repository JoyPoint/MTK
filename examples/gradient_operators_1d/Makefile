# Makefile for example: gradient_operators_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

gradient_operators_1d.exe: \
  gradient_operators_1d.cc \
  gradient_operators_1d.o
	$(LINKER) gradient_operators_1d.o $(LIBS) \
  -o gradient_operators_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

gradient_operators_1d.o: gradient_operators_1d.cc $(MTK_LIB)
	$(CC) gradient_operators_1d.cc $(CCFLAGS) -o gradient_operators_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) gradient_operators_1d.exe >& \
          gradient_operators_1d.memcheck

clean:
	rm -f *~ *.o *.tex *.exe
