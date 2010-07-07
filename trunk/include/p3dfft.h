/*
! This file is part of P3DFFT library
!
! Version 2.3
!
! Copyright (C) 2006-2008 Dmitry Pekurovsky
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!----------------------------------------------------------------------------
*/

#include <stdlib.h>

#ifdef IBM

#define FORT_MOD_NAME(NAME) __p3dfft_NMOD_##NAME
#define FORTNAME(NAME) NAME

#elif defined PGI

#define FORT_MOD_NAME(NAME) p3dfft_##NAME##_
#define FORTNAME(NAME) NAME##_

#elif defined GCC

#define FORT_MOD_NAME(NAME) __p3dfft__##NAME
#define FORTNAME(NAME) NAME##_

#else

#define FORT_MOD_NAME(NAME) p3dfft_mp_##NAME##_
#define FORTNAME(NAME) NAME##_

#endif

extern void FORT_MOD_NAME(p3dfft_setup)(int *dims,int *nx,int *ny,int *nz, int *ow);
extern void FORT_MOD_NAME(p3dfft_get_dims)(int *,int *,int *,int *);

#ifndef SINGLE_PREC
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(double *A,double *B);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(double *A,double *B);
#else
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(float *A,float *B);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(float *A,float *B);
#endif

extern void FORT_MOD_NAME(p3dfft_clean)();
extern void FORTNAME(abort)();

extern void p3dfft_setup(int *dims,int nx,int ny,int nz,int ovewrite);
extern void p3dfft_get_dims(int *,int *,int *,int );

#ifndef SINGLE_PREC
extern void p3dfft_ftran_r2c(double *A,double *B);
extern void p3dfft_btran_c2r(double *A,double *B);
#else
extern void p3dfft_ftran_r2c(float *A,float *B);
extern void p3dfft_btran_c2r(float *A,float *B);
#endif

extern void p3dfft_clean();

inline void p3dfft_setup(int *dims,int nx,int ny,int nz, int overwrite)
{
  FORT_MOD_NAME(p3dfft_setup)(dims,&nx,&ny,&nz,&overwrite);
}

inline void p3dfft_get_dims(int *start,int *end,int *size,int conf)
{
  FORT_MOD_NAME(p3dfft_get_dims)(start,end,size,&conf);
}

#ifndef SINGLE_PREC
inline void p3dfft_ftran_r2c(double *A,double *B)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c)(A,B);
}
#else
inline void p3dfft_ftran_r2c(float *A,float *B)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c)(A,B);
}

#endif

#ifndef SINGLE_PREC
inline void p3dfft_btran_c2r(double *A,double *B)
{
  FORT_MOD_NAME(p3dfft_btran_c2r)(A,B);
}
#else
inline void p3dfft_btran_c2r(float *A,float *B)
{
  FORT_MOD_NAME(p3dfft_btran_c2r)(A,B);
}
#endif

inline void p3dfft_clean()
{
  FORT_MOD_NAME(p3dfft_clean)();
}

inline void FORTNAME(abort)() 
{
  abort();
}
/* config.h.  Generated by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* arguments passed to configure script */
#define CONFIGURE_ARGS "'--enable-pgi' '--enable-fftw' '--with-fftw=/opt/fftw/3.2.1/pgi' 'FCFLAGS=-fastsse -tp barcelona-64 -Mextend -byteswapio' 'CFLAGS=-fastsse -tp barcelona-64' 'LDFLAGS=-lmpi_f90 -lmpi_f77 -lmyriexpress'"

/* Define if you want to enable C convention for processor dimensions */
/* #undef DIMS_C */

/* Define if you want to use the ESSL library instead of FFTW */
/* #undef ESSL */

/* Define whether you want to enable estimation */
/* #undef ESTIMATE */

/* Define if you want to use the FFTW library */
#define FFTW 1

/* Define if you want to compile P3DFFT using GCC compiler */
/* #undef GCC */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if you want to compile P3DFFT using IBM compiler */
/* #undef IBM */

/* Define if you want to compile P3DFFT using Intel compiler */
/* #undef INTEL */

/* Define if you want to enable the measure algorithm */
#define MEASURE 1

/* Define if you want to override the default value of NBL_X */
/* #undef NBL_X */

/* Define if you want to override the default value of NBL_Y */
/* #undef NBL_Y */

/* Define if you want 1D decomposition */
/* #undef ONED */

/* Name of package */
#define PACKAGE "p3dfft"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "dmitry@sdsc.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "P3DFFT"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "P3DFFT 2.3.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "p3dfft"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.3.1"

/* Define if you want to enable the patient algorithm */
/* #undef PATIENT */

/* Define if you want to compile P3DFFT using PGI compiler */
#define PGI 1

/* Define if you want to compile P3DFFT in single precision */
/* #undef SINGLE_PREC */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if you want to enable stride-1 data structures */
/* #undef STRIDE1 */

/* Define if you want to MPI_Alltoall instead of MPI_Alltotallv */
/* #undef USE_EVEN */

/* Version number of package */
#define VERSION "2.3.1"

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */
