FC = mpif90
AM_FCFLAGS = $(FCFLAGS) $(DFLAGS) $(INCLUDES) 
INCLUDES = -I$(top_builddir)/include/
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW3_LIB)/libfftw3.a

noinst_PROGRAMS =  test_sine.x test_sine_inplace.x test_rand.x test_spec.x

test_sine_x_SOURCES = driver_sine.f90

test_sine_inplace_x_SOURCES =  driver_sine_inplace.f90

test_rand_x_SOURCES = driver_rand.f90

test_spec_x_SOURCES = driver_spec.f90

clean-local:
	-test -z "*.x" || rm -f *.x
