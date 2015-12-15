# Makefile for example: difussion_3d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

difussion_3d.exe: difussion_3d.cc difussion_3d.o
	$(LINKER) difussion_3d.o $(LIBS) -o difussion_3d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

difussion_3d.o: difussion_3d.cc $(MTK_LIB)
	$(CC) difussion_3d.cc $(CCFLAGS) -o difussion_3d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) difussion_3d.exe >& difussion_3d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
