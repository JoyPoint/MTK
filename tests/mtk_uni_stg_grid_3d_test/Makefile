# Makefile for test: mtk_uni_stg_grid_3d_test.cc

include ../../Makefile.conf

NAME = mtk_uni_stg_grid_3d_test

#   Linking statements.
#   __________________________________________________________________

$(NAME).exe: $(NAME).cc $(NAME).o
	$(LINKER) $(NAME).o $(LIBS) -o $(NAME).exe
	@echo ' '

#   Compiling statements.
#   __________________________________________________________________

$(NAME).o: $(NAME).cc $(MTK_LIB)
	$(CC) $(NAME).cc $(CCFLAGS) -o $(NAME).o

#   Utilities.
#   __________________________________________________________________

memcheck:
	valgrind $(MEMCHECK_OPTS) $(NAME).exe >& $(NAME).memcheck

clean:
	rm -f *~ *.o *.exe *.dat *.png *.tex *.pdf
