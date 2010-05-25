
!========================================================
!  3D FFT Forward transform with 2D domain decomposition
!
!  The order of array elements in memory is the same in all stages: (x,y,z) 
!
! Input: XgYZ - an array of real, x dimension is entirely local within 
! processors memory  while y and z are block-distributed among processors 
! in 2D grid
! Output: XYZg - comlpex array, with entire Z dimension local,
!               while x and y are block-distributed among processors in 2D grid
!
! The arrays may occupy the same memory space
! In this case their first elements should coincide. 
! Naturally, output overwrites input
!
      subroutine p3dfft_ftran_r2c (XgYZ,XYZg)
!========================================================

      use fft_spec
      implicit none

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif
      integer x,y,z,i,err,nx,ny,nz
      integer(8) Nl
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft
!
! FFT transform (R2C) in X for all z and y
!
      if(jisize * kjsize .gt. 0) then
         call init_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize) 

         timers(5) = timers(5) - MPI_Wtime()
         call exec_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)
         timers(5) = timers(5) + MPI_Wtime()
      endif

! Exchange data in rows 
      if(iproc .gt. 1) call fcomm1(buf,buf,timers(1))

! FFT transform (C2C) in Y for all x and z, one Z plane at a time

      if(iisize * kjsize .gt. 0) then
         call init_f_c(buf,iisize,1,buf,iisize,1,ny,iisize)
         timers(6) = timers(6) - MPI_Wtime()
         do z=1,kjsize
            call ftran_y_zplane(buf,z-1,iisize,kjsize,iisize,1, &
                 buf,z-1,iisize,kjsize,iisize,1,ny,iisize)
         enddo
         timers(6) = timers(6) + MPI_Wtime()

      endif

! Exchange data in columns
      if(jproc .gt. 1) then
#ifdef STRIDE1
! For stride1 option combine second transpose with transform in Z
         call init_f_c(buf,1,nz, buf,1,nz,nz,iisize*NB)
         call fcomm2_trans(buf,XYZg,timers(2),timers(7))
#else


! FFT Transform (C2C) in Z for all x and y

! Transpose y-z

         call fcomm2(buf,XYZg,timers(2),timers(7))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

         if(iisize * jjsize .gt. 0) then
            call init_f_c(XYZg,iisize*jjsize, 1, &
                 XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         
            timers(7) = timers(7) - MPI_Wtime()
            call exec_f_c2(XYZg,iisize*jjsize, 1, &
                 XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
            timers(7) = timers(7) + MPI_Wtime()
         endif

#endif

      else
#ifdef STRIDE1

         call reorder_trans_f(buf,XYZg)
#else
         Nl = iisize*jjsize*nz
         call ar_copy(buf,XYZg,Nl)
         call init_f_c(XYZg,iisize*jjsize, 1, &
              XYZg,iisize*jjsize, 1,nz,iisize*jjsize)         
         call exec_f_c2(XYZg,iisize*jjsize, 1, &
              XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
#endif
      endif


      return
      end subroutine

