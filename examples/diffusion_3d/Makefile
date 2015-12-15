# Makefile for example: diffusion_3d.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

diffusion_3d.exe: diffusion_3d.cc diffusion_3d.o
	$(LINKER) diffusion_3d.o $(LIBS) -o diffusion_3d.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

diffusion_3d.o: diffusion_3d.cc $(MTK_LIB)
	$(CC) diffusion_3d.cc $(CCFLAGS) -o diffusion_3d.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) diffusion_3d.exe >& diffusion_3d.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
