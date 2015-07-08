## Installing P3DFFT ##
The latest version of P3DFFT can be found here. Once you have extracted the package, you must take the following steps to complete the setup:

  1. Run the configure script
  1. Run 'make'
  1. Run 'make install'

# How to compile P3DFFT #
P3DFFT uses a configure script to create Makefiles for compiling the library as well as several examples. This page will step you through the process of running the configure script properly.

Run "configure --help" for complete list of options.
Recommended options: --enable-stride1. For Cray XT platforms also recommended --enable-useeven.

Currently the package supports four compiler suites: PGI, Intel, IBM and GNU. Some examples of compiling on several systems are given below. Users may need to customize
as needed. If you wish to share more examples or to request or contribute in support for other compilers, please write to dmitry@sdsc.edu. If you give us an account on your system we will work with you to customize the installation.

| **Argument** | **Notes** | **Description** | **Example** |
|:-------------|:----------|:----------------|:------------|
| --prefix=PREFIX | Mandatory for users without access to /usr/local | This argument will install p3dfft to PREFIX when you run make install. By default, configure will install to /usr/local | --prefix=$HOME/local/ |
| --enable-gnu, --enable-ibm, --enable-intel, --enable-pgi, --enable-cray | Mandatory | These arguments will prepare p3dfft to be built by a specific compiler. You must only choose one option. | --enable-pgi |
| --enable-fftw, --enable-essl | Mandatory | These arguments will prepare p3dfft to be used with either the FFTW or ESSL library. You must only choose one option. | --enable-fftw |
| --with-fftw= _FFTWLOCATION_ | Mandatory if --enable-fftw is used | This argument specifies the path location for the FFTW library; it is mandatory if you are planning to use p3dfft with the FFTW library. | --enable-fftw --with-fftw=$FFTW\_HOME |
| --enable-oned | Optional, | This argument is for 1D decomposition. The default is 2D decomposition but can be made to 1D by setting up a grid 1xn when running the code. | --enable-oned |
| --enable-estimate | Optional, use only with --enable-fftw | If this argument is passed, the FFTW library will not use run-time tuning to select the fastest algorithm for computing FFTs. | --enable-estimate |
| --enable-measure | Optional, enabled by default, use only with --enable-fftw | For search-once-for-the-fast algorithm (takes more time on p3dfft\_setup()). | --enable-measure |
| --enable-patient | Optional, use only with --enable-fftw | For search-once-for-the-fastest-algorithm (takes much more time on p3dfft\_setup()). | --enable-patient |
| --enable-dimsc | Optional  | To assign processor rows and columns in the Cartesian processor grid according to C convention. The default is Fortran convention which is recommended. This option does not affect the order of storage of arrays in memory. | --enable-dimsc |
| --enable-useeven | Optional, recommended for Cray XT | This argument is for using MPI\_Alltoall instead of MPI\_Alltotallv. This will pad the send buffers with zeros to make them of equal size; not needed on most architecture but may lead to better results on Cray XT. | --enable-useeven |
| --enable-stride1 | Optional, recommended | To enable stride-1 data structures on output (this may in some cases give some advantage in performance). You can define loop blocking factors NLBX and NBLY to experiment, otherwise they are set to default values. | --enable-stride1 |
| --enable-nblx | Optional  | To define loop blocking factor NBL\_X | --enable-nblx=32 |
| --enable-nbly1 | Optional  | To define loop blocking factor NBL\_Y1 | --enable-nbly1=32 |
| --enable-nbly2 | Optional  | To define loop blocking factor NBL\_Y2 | --enable-nbly2=32 |
| --enable-nblz | Optional  | To define loop blocking factor NBL\_Z | --enable-nblz=32 |
| --enable-single | Optional  | This argument will compile p3dfft in single-precision. By default, configure will setup p3dfft to be compiled in double-precision. | --enable-single |
| FC=<Fortran compiler> | Optional  | Fortran compiler | FC=mpfort   |
| FCFLAGS="<Fortran compiler flags>" | Optional, recommended | Fortran compiler flags | FCFLAGS="-Mextend" |
| CC=<C compiler> | Optional  | C compiler      | CC=mpcc     |
| CFLAGS="<C compiler flags>" | Optional, recommended | C compiler flags | CFLAGS="-fastsse" |
| LDFLAGS="<linker flags>" | Mandatory (depending on platform) | Linker flags    | LDFLAGS=-lmpi\_f90 -lmpi\_f77" |

# Compiling on Gordon/Trestles (SDSC) #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| PGI          | pgi, fftw, mvapich2\_ib | ./configure --enable-pgi --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc |
| Intel        | intel, fftw, mvapich2\_ib | ./configure --enable-intel --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc|
| GNU          | gnu, fftw, mvapich2\_ib | ./configure --enable-intel --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc |


# Compiling on Triton (SDSC) #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| PGI          | pgi, fftw, mpi | ./configure --enable-pgi --enable-fftw --with-fftw="$FFTWHOME" FCFLAGS="-fastsse -tp barcelona-64 -Mpreprocess" CFLAGS="-fastsse -tp barcelona-64" LDFLAGS="-lmpi\_f90 -lmpi\_f77 -lmyriexpress" |
| Intel        | intel, fftw, mpi | ./configure --enable-intel --enable-fftw --with-fftw="$FFTWHOME" FCFLAGS="-O3 -xW -132 -fpp" CFLAGS="-O3 -xW" LDFLAGS="-lmpi\_f90 -lmpi\_f77 -lmyriexpress -limf" |
| GNU          | gnu, fftw, mpi | ./configure --enable-gnu --enable-fftw --with-fftw="$FFTWHOME" LDFLAGS="-lmpi\_f90 -lmpi\_f77" |

# Compiling on IBM Power7 #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| mpfort       | essl        | ./configure FC=mpfort FCFLAGS="-qcclines -qarch=pwr7 -qstrict -qnosave -qtune=pwr7 -qhot -qsimd=auto -qcache=auto -qsmp=omp -qthreaded -O3 -q64 -qfloat=hsflt:fltint" CC=mpcc CFLAGS="-q64 -DNUS\_XCOMP" CPP=/usr/bin/cpp --enable-essl --enable-ibm --host="bd-login" |

# Compiling on Stampede (TACC) #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| Intel        | fftw3       | ./configure --enable-fftw --enable-intel --with-fftw=$TACC\_FFTW3\_DIR FC=mpif90 CC=mpicc|

# Compiling on Hopper (NERSC) #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| PGI          | PrgEnv-pgi, fftw | ./configure --enable-pgi --enable-fftw --with-fftw=$FFTW\_DIR/.. FC=ftn CC=cc |
| Intel        | PrgEnv-intel, fftw| ./configure --enable-intel --enable-fftw --with-fftw=$FFTW\_DIR/.. FC=ftn CC=cc|
| GNU          | PrgEnv-gnu, fftw | ./configure --enable-gnu --enable-fftw --with-fftw=$FFTW\_DIR/.. FC=ftn CC=cc|
| Cray         | PrgEnv-cray, fftw | ./configure --enable-cray --enable-fftw --with-fftw=$FFTW\_DIR/.. FC=ftn CC=cc|

# Compiling on Blacklight (PSC) #
| **Compiler** | **Modules** | **Arguments** |
|:-------------|:------------|:--------------|
| Intel        | fftw        | ./configure --enable-fftw --enable-intel --with-fftw=$FFTW\_LIB/.. CC=icc FC=ifort FCFLAGS=-lmpi|

# Compiling on Mira/Cetus/Vesta (ALCF) #
| **Compiler** | **Arguments** |
|:-------------|:--------------|
| IBM XL       | ./configure --enable-ibm --enable-essl --with-essl=/soft/libraries/essl/current FC=mpixlf90\_r CC=mpixlc\_r|
| GNU          | ./configure --enable-gnu --enable-fftw --with-fftw=/soft/libraries/alcf/current/{xl,gcc}/FFTW3 FC=mpif90 CC=mpicc|