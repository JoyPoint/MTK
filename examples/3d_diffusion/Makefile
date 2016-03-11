# Makefile for example: 3d_diffusion.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

3d_diffusion.exe: 3d_diffusion.cc 3d_diffusion.o
	$(LINKER) 3d_diffusion.o $(LIBS) -o 3d_diffusion.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

3d_diffusion.o: 3d_diffusion.cc $(MTK_LIB)
	$(CC) 3d_diffusion.cc $(CCFLAGS) -o 3d_diffusion.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 3d_diffusion.exe >& 3d_diffusion.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
