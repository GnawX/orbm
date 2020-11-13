# Top level makefile for qe-orbm

all:    build

build:
	@echo "Building qe-orbm..."
	$(MAKE) -C src

clean:
	@echo "Cleaning qe-orbm..."
	if test -s src/Makefile ; then ( $(MAKE) -C src clean ); fi
	-/bin/rm -f bin/orbm.x

distclean:
	$(MAKE) -C src distclean
	-/bin/rm -f config.log config.status makedeps.sh
	-/bin/rm -Rf autom4te.cache

