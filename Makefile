## Main makefile.

include Makefile.inc

all: mtklib test example

mtklib:
	@( cd $(SRC); $(MAKE) )

test:
	@( cd $(TESTS); $(MAKE) )

example:
	@( cd $(EXAMPLES); $(MAKE) )

gendoc:
	@( $(DOCGEN) $(DOCFILE) )
	@( cd $(DOC)/latex; $(MAKE) )

clean: cleanroot cleanlib cleantest cleanexample

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

bye: clean gendoc

help:
	@echo '-----'
	@echo 'Makefile for the MTK.'
	@echo ' '
	@echo 'Options are:'
	@echo '- all: builds he library, the tests, and examples.'
	@echo '- mtklib: builds the library.'
	@echo '- test: builds the test files.'
	@echo '- example: builds the examples.'
	@echo ' '
	@echo '- gendoc: generates the documentation for the library.'
	@echo ' '
	@echo '- clean: cleans ALL the generated files.'
	@echo '- cleanlib: cleans the generated archive and object files.'
	@echo '- cleantest: cleans the generated tests executables.'
	@echo '- cleanexample: cleans the generated examples executables.'
	@echo '-----'
