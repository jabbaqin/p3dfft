F77 = $(FC)
FC = mpif90
AM_FFLAGS = $(FCFLAGS) $(DFLAGS) #-fastsse -tp barcelona-64 -Mextend -byteswapio
INCLUDES = -I$(FFTW3_INC)

lib_LIBRARIES = libp3dfft.a

libp3dfft_a_SOURCES = fft_spec.F module.F fft_init.F fft_exec.F wrap.f

# all-local:
#	-mv libp3dfft.a $(top_builddir)/lib/
#	-mv p3dfft.mod $(top_builddir)/include/

clean-local:
	-test -z "*.mod" || rm -f *.mod
