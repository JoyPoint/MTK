## Main makefile.

include Makefile.conf

all: mtklib test example

mtklib:
	@( cd $(SRC); $(MAKE) )

test:
	@( cd $(TESTS); $(MAKE) )

example:
	@( cd $(EXAMPLES); $(MAKE) )

testall:
	@( cd $(TESTS); $(MAKE) testall )

gendoc:
	@( $(DOCGEN) $(DOCFILE) )
	@( cd $(DOC)/latex; $(MAKE) )

clean: cleanroot cleanlib cleantest # cleanexample

cleanroot:
	@(rm -f *~)

cleanlib:
	@( cd $(SRC); $(MAKE) clean )
	@( cd $(LIB); $(MAKE) clean )

cleantest:
	@( cd $(TESTS); $(MAKE) clean )

cleanexample:
	@( cd $(EXAMPLES); $(MAKE) clean )

memcheck:
	@( cd $(TESTS); $(MAKE) memcheck)
	@( cd $(EXAMPLES); $(MAKE) memcheck)

help:
	@echo '-----'
	@echo 'Makefile for the MTK.'
	@echo ' '
	@echo 'Options are:'
	@echo '- all: builds the library, the tests, and examples.'
	@echo '- mtklib: builds the library.'
	@echo '- test: builds the test files.'
	@echo '- example: builds the examples.'
	@echo ' '
	@echo '- testall: runs all the tests.'
	@echo ' '
	@echo '- gendoc: generates the documentation for the library.'
	@echo ' '
	@echo '- clean: cleans all the generated files.'
	@echo '- cleanlib: cleans the generated archive and object files.'
	@echo '- cleantest: cleans the generated tests executables.'
	@echo '- cleanexample: cleans the generated examples executables.'
	@echo '-----'
