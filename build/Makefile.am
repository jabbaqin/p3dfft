F77 = $(FC)
FC = mpif90
AM_FFLAGS = $(FCFLAGS) $(DFLAGS) 
INCLUDES = -I$(FFTW3_INC)

lib_LIBRARIES = libp3dfft.a

libp3dfft_a_SOURCES = fft_spec.F90 module.F90 fft_init.F90 fft_exec.F90 wrap.f90

clean-local:
	-test -z "*.mod" || rm -f *.mod
