# Makefile for example: mimetic_threshold_div_grad.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

mimetic_threshold_div_grad.exe: \
  mimetic_threshold_div_grad.cc mimetic_threshold_div_grad.o
	$(LINKER) mimetic_threshold_div_grad.o $(LIBS) \
  -o mimetic_threshold_div_grad.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

mimetic_threshold_div_grad.o: mimetic_threshold_div_grad.cc $(MTK_LIB)
	$(CC) mimetic_threshold_div_grad.cc $(CCFLAGS) \
  -o mimetic_threshold_div_grad.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) mimetic_threshold_div_grad.exe >& \
  mimetic_threshold_div_grad.memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.eps *.tex *.aux *.log *.pdf
