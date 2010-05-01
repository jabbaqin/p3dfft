SHELL=/bin/csh

include makeinc

export P3DFFT_ROOT = ../..
export P3DFFT_LIB = p3dfft
export P3DFFT_LOC = $(P3DFFT_ROOT)/lib
export P3DFFT_INC = $(P3DFFT_ROOT)/include

.PHONY: C FORTRAN clean install

all:
	cd FORTRAN && make -f makefile.loc
	cd C && make -f makefile.loc
C:
	cd C && make -f makefile.loc
FORTRAN:
	cd FORTRAN && make -f makefile.loc
clean:
	cd FORTRAN && make clean -f makefile.loc
	cd C && make clean -f makefile.loc
install:
	cp -f C/test_rand.x test_rand.C.x
	cp -f C/test_sine.x test_sine.C.x
	cp -f C/test_sine_inplace.x test_since_inplace.C.x
	cp -f C/test_spec.x test_spec.C.x
	cp -f FORTRAN/test_rand.x test_rand.F.x
	cp -f FORTRAN/test_sine.x test_sine.F.x
	cp -f FORTRAN/test_sine_inplace.x test_sine_inplace.F.x
	cp -f FORTRAN/test_spec.x test_spec.F.x
