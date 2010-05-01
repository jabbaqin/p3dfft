CC = mpicc
AM_CPPFLAGS = $(DFLAGS) #-fastsse -tp barcelona-64 -DPGI 
#AM_LDFLAGS = -lmpi_f90 -lmpi_f77 -pgf90libs -lmyriexpress 
INCLUDES = -I$(top_builddir)/include/ 
LDADD = $(top_builddir)/lib/libp3dfft.a $(FFTW3_LIB)/libfftw3.a

noinst_PROGRAMS =  test_sine.x test_sine_inplace.x test_rand.x test_spec.x

test_sine_x_SOURCES = driver_sine.c

test_sine_inplace_x_SOURCES =  driver_sine_inplace.c

test_rand_x_SOURCES = driver_rand.c

test_spec_x_SOURCES = driver_spec.c

clean-local:
	-test -z "*.x" || rm -f *.x
