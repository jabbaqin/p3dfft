#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(P3DFFT, 2.3.1, dmitry@sdsc.edu)
AM_INIT_AUTOMAKE(no-dependencies)
AC_CONFIG_SRCDIR([build/module.F])
AC_CONFIG_HEADER([config.h])

# check whether to enable single precision
AC_MSG_CHECKING([whether to enable single precision])
AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single], [compile P3DFFT in single precision (default is double precision)])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
	AC_SUBST(DFLAGS, "$DFLAGS -DSINGLE_PREC")
fi

# check whether to enable 1D decomposition
AC_MSG_CHECKING([whether to enable 1D decomposition])
AC_ARG_ENABLE(oned, [AC_HELP_STRING([--enable-oned], [for 1D decomposition (default is 2D but it can be made 1D through setting up grid 1 x n when running the code) ])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DONED")

fi

# check whether to use ESSL
AC_MSG_CHECKING([whether to use ESSL library])
AC_ARG_ENABLE(essl, [AC_HELP_STRING([--enable-essl], [for using ESSL (on IBM) in 1D FFTs.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "DFLAGS -DESSL")
fi

# check whether to enable estimation
AC_MSG_CHECKING([whether to enable estimation])
AC_ARG_ENABLE(estimate, [AC_HELP_STRING([--enable-estimate], [for initialization thoroughness in setting up FFTW calls (this option only goes with --enable-fftw) .])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DESTIMATE")
fi

# check whether to enable the measure algorithm
AC_MSG_CHECKING([whether to enable the emasure algorithm])
AC_ARG_ENABLE(measure, [AC_HELP_STRING([--enable-measure], [for search-once-for-the-fast-algorithm (takes more time on p3dfft_setup()).])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DMEASURE")
fi

# check whether to enable the patient algorithm
AC_MSG_CHECKING([whehter to enable the patient algorithm])
AC_ARG_ENABLE(patient, [AC_HELP_STRING([--enable-patient], [for search-once-for-the-fastest-algorithm (takes much more time on p3dfft_setup()).])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DPATIENT")
fi

# check whether to enable C convention for processor dimensions
AC_MSG_CHECKING([whether to enable C convention for processor dimensions])
AC_ARG_ENABLE(dimsc, [AC_HELP_STRING([--enable-dimsc], [to assign processor rows and columns according to C convention, rather than "Fortran", which is default.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DDIMS_C")
fi

# check whether to reverse processor dimensions
AC_MSG_CHECKING([whether to reverse processor dimensions])
AC_ARG_ENABLE(reversedims, [AC_HELP_STRING([--enable-reversedims], [to interchange processor rows and columns from default.])], ok=$enableval, ok=yes)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
	AC_SUBST(DFLAGS, "$DFLAGS -DREVERSE_DIMS")
fi

# check whether to use MPI_Alltoall instead of MPI_Alltotallv
AC_MSG_CHECKING([whether to use MPI_Alltoall instead of MPI_Alltotallv])
AC_ARG_ENABLE(useeven, [AC_HELP_STRING([--enable-useeven], [for using MPI_Alltoall instead of MPI_Alltotallv. This method pads the send buffers with zeros to make them equal size. This options is not needed on most architectures.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DUSE_EVEN")
fi

# check whether to enable stride-1 data structures
AC_MSG_CHECKING([whether to enable stride-1 data structures])
AC_ARG_ENABLE(stride1, [AC_HELP_STRING([--enable-stride1], [to enable stride-1 data structures on output (this may in some cases give some advantage in performance). You can define loop blocking factors NBL_X and NBL_Y to experiment, otherwise they are set to default values.])], ok=$enableval, ok=no)
AC_MSG_RESULT([$ok])
if test "$ok" = "yes"; then
        AC_SUBST(DFLAGS, "$DFLAGS -DSTRIDE1")
fi

#check for fftw3
AC_ARG_WITH(fftw3, [AC_HELP_STRING([--with-fftw3], [path to location of FFTW3 library])], loc=$withval, loc="")
if test "$loc" != ""; then
        AC_SUBST(DFLAGS, "$DFLAGS -DFFTW")
	AC_SUBST(FFTW3_INC, "$loc/include")
	AC_SUBST(FFTW3_LIB, "$loc/lib")
	
fi


# Checks for programs.
AC_PROG_FC(mpif90)
AC_CHECK_PROGS(CC, mpicc, $CC)
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile
                 build/Makefile
		 sample/Makefile
		 sample/C/Makefile
		 sample/FORTRAN/Makefile])

AC_SUBST(TESTLIBS, $TESTLIBS)
AC_OUTPUT
