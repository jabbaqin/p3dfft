F77 = $(FC)
FC = mpif90
AM_FFLAGS = $(FCFLAGS) $(DFLAGS) $(INCLUDES) 
INCLUDES = -I$(top_builddir)/include/
LDADD = $(top_builddir)/build/libp3dfft.a $(FFTW3_LIB)/libfftw3.a

noinst_PROGRAMS =  test_sine.x test_sine_inplace.x test_rand.x test_spec.x

test_sine_x_SOURCES = driver_sine.f

test_sine_inplace_x_SOURCES =  driver_sine_inplace.f

test_rand_x_SOURCES = driver_rand.f

test_spec_x_SOURCES = driver_spec.f

clean-local:
	-test -z "*.x" || rm -f *.x
