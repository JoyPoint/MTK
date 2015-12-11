# Makefile for example: curl_2d_angular_velocity.exe

include ../../Makefile.inc

#   Linking statements.
#   __________________________________________________________________

curl_2d_angular_velocity.exe: curl_2d_angular_velocity.cc \
  curl_2d_angular_velocity.o
	$(LINKER) curl_2d_angular_velocity.o $(LIBS) \
          -o curl_2d_angular_velocity.exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

curl_2d_angular_velocity.o: curl_2d_angular_velocity.cc $(MTK_LIB)
	$(CC) curl_2d_angular_velocity.cc $(CCFLAGS) \
         -o curl_2d_angular_velocity.o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) curl_2d_angular_velocity.exe >& \
          curl_2d_angular_velocity.memcheck

clean:
	rm -f *~ *.o *.dat *.png *.exe
