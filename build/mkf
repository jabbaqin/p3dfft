# This file is part of P3DFFT library
#
# Version 2.3
#
# Copyright (C) 2006-2008 Dmitry Pekurovsky
#
#    P3DFFT is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#----------------------------------------------------------------------------

SHELL=/bin/csh

include ./makeinc

P3DFFT_ROOT = ..
P3DFFT_LIB = libp3dfft.a

# Be sure to set DFLAGS options in makeinc. The menu is as follows

# -DDOUBLE_PREC : for double precision (default is double precision)
#
# -DONED : for 1D decomposition (default is 2D but it can be made 1D through
#    setting up grid 1 x n when running the code) 
#
# -DFFTW or -DESSL for using FFTW3 library or ESSL (on IBM) in 1D FFTs.
#   One of these has to be set.
#
# -DESTIMATE, -DMEASURE or -DPATIENT (default -DMEASURE) for initialization thouroughness in setting up FFTW calls   (this option only goes with -DFFTW)
# For search-once-for-a-fast-algorithm "-DMEASURE" (takes more time 
#on p3dfft_setup())
# For search-once-for-the-fastest-algorithm "-DPATIENT" (takes much 
#more time on p3dfft_setup())
#
# -DREVERSE_DIMS to interchange processor rows and columns from default 
# -DUSE_EVEN for using MPI_Alltoall instead of MPI_Alltoallv. This method
#   pads the send buffers with zeros to make them equal size. This options
#   is not needed on most architectures.


# ----------------

FFT3DLIB = module.o fft_spec.o fft_init.o fft_exec.o wrap.o

all: lib install
lib: $(FFT3DLIB)
	$(AR) $(ARFLAGS) $(P3DFFT_LIB) $(FFT3DLIB)	
install: 
	if(!(-e $(P3DFFT_ROOT))) mkdir $(P3DFFT_ROOT)
	if (!(-e $(P3DFFT_ROOT)/lib)) mkdir $(P3DFFT_ROOT)/lib	
	cp $(P3DFFT_LIB) $(P3DFFT_ROOT)/lib
	if (!(-e $(P3DFFT_ROOT)/include)) mkdir $(P3DFFT_ROOT)/include
	cp p3dfft.mod $(P3DFFT_ROOT)/include
	if (!(-e $(P3DFFT_ROOT)/include/p3dfft.h)) -cp ../include/p3dfft.h $(P3DFFT_ROOT)/include

module.o: module.F fft_spec.o
fft_init.o: fft_init.F module.o fft_spec.o
fft_exec.o: fft_exec.F module.o fft_spec.o
wrap.o: wrap.f module.o

.c.o: 
		$(CC) -c $(CFLAGS) *.c
.F.o:
		$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o: 
		$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm $(FFT3DLIB) *.mod




