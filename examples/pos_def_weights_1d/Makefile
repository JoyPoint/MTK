# Makefile for example: pos_def_weights_1d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

pos_def_weights_1d.exe: pos_def_weights_1d.cc pos_def_weights_1d.o
	$(LINKER) pos_def_weights_1d.o $(LIBS) -o pos_def_weights_1d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

pos_def_weights_1d.o: pos_def_weights_1d.cc $(MTK_LIB)
	$(CC) pos_def_weights_1d.cc $(CCFLAGS) -o pos_def_weights_1d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) pos_def_weights_1d.exe >& \
          pos_def_weights_1d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
