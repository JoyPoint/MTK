# Makefile for example: 2d_angular_velocity.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

2d_angular_velocity.exe: 2d_angular_velocity.cc 2d_angular_velocity.o
	$(LINKER) 2d_angular_velocity.o $(LIBS) -o 2d_angular_velocity.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

2d_angular_velocity.o: 2d_angular_velocity.cc $(MTK_LIB)
	$(CC) 2d_angular_velocity.cc $(CCFLAGS) -o 2d_angular_velocity.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) 2d_angular_velocity.exe >& \
          2d_angular_velocity.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
