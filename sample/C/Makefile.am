CC = mpicc
AM_CPPFLAGS = $(DFLAGS) 
INCLUDES = -I$(top_builddir)/include/ 
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW3_LIB)/libfftw3.a $(FFTW3_LIB)/libfftw3f.a

noinst_PROGRAMS =  test_sine.x test_sine_inplace.x test_rand.x test_spec.x

test_sine_x_SOURCES = driver_sine.c

test_sine_inplace_x_SOURCES =  driver_sine_inplace.c

test_rand_x_SOURCES = driver_rand.c

test_spec_x_SOURCES = driver_spec.c

clean-local:
	-test -z "*.x" || rm -f *.x
