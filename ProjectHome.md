## User Survey is now open (until September 15, 2014) ##

Please fill out the [User Survey](https://docs.google.com/forms/d/1mOVt1QHEe16rj6jHH-0TAw9wLg9QfMnMprASYAgapCY/viewform) to provide feedback to the authors and
request features.

## P3DFFT introduction ##
Parallel Three-Dimensional Fast Fourier Transforms, dubbed P3DFFT, is a library for large-scale computer simulations on parallel platforms. 3D FFT is an important algorithm for simulations in a wide range of fields, including studies of turbulence, climatology, astrophysics and material science. This project was initiated at [San Diego Supercomputer Center (SDSC)](http://www.sdsc.edu) at UC San Diego by its main author Dmitry Pekurovsky, Ph.D.

P3DFFT uses 2D, or pencil, decomposition. This overcomes an important limitation to scalability inherent in FFT libraries implementing 1D (or slab) decomposition: the number of processors/tasks used to run this problem in parallel can be as large as N^2, were N is the linear problem size. This approach has shown good scalability up to 131,072 cores. P3DFFT is optimized for large data sets.

P3DFFT is written in Fortran90/MPI and is optimized for parallel performance. C interface is available, as are detailed documentation and examples in both Fortran and C.
A configure script is supplied for ease of installation. This package depends on a serial
FFT library such as [Fastest Fourier Transform in the West](http://fftw.org/) (FFTW) or IBM's [ESSL](http://www-03.ibm.com/systems/software/essl/index.html).

In the forward transform, given an input of an array of 3D real values, an output of 3D complex array of Fourier coefficients is returned. Current features include:
  1. real-to-complex/complex-to-real FFT in 3D
  1. real-to-complex FFT in 2D followed by sine/cosine/Chebyshev/empty transform, and the reverse for backward transform.
  1. pruned transforms (less than full input or output)
  1. in-place or out-of-place transforms
  1. multi-variable transforms
For more information see the [User Guide](UserGuide.md). For installation instructions [click here](Install.md).

## Citation information ##

Please acknowledge/cite use of P3DFFT as follows:
D. Pekurovsky, “P3DFFT: a framework for parallel computations of Fourier transforms in three dimensions”, SIAM Journal on Scientific Computing 2012, Vol. 34, No. 4, pp. C192-C209. This paper can be obtained [here](http://p3dfft.googlecode.com/svn/docs/P3DFFT_SISC_final.pdf).

## Version History ##

2.7.1 - Added multi-variable transforms (p3dfft\_ftran\_r2c\_many, p3dfft\_btran\_c2r\_many)

2.6.1 - Added pruned transforms

> Added user-defined communicator

2.5.1 - Added cosine/sine/Chebyshev/empty transform in addition to Fourier

## Communication ##

Be sure to subscribe to the [project mailing list](http://groups.google.com/group/p3dfft) where you can discuss topics of interest with other users and developers and get timely news regarding this library. Click on Issues tab above to see the list of outstanding problems or report a new problem or suggestion. You can also reach the main author Dmitry Pekurovsky at dmitry@sdsc.edu.

We are interested to hear your experiences, and suggestions for future releases. Please let us know if you are interested in contributing to future development of the library.
The author is also interested in collaborating on computational science projects where 3D FFT and related transforms are applied.

## Future work ##

  * Hybrid MPI/OpenMP
  * One-sided communication
  * Generalized data layout options