!========================================================
!  3D FFT inverse transform with 2D domain decomposition
!
!  The order of array elements in memory is the same in all stages: (x,y,z) 
!
! Input: XYZg - comlpex array, with entire Z dimension local,
!               while x and y are block-distributed among processors in 2D grid
! Output: XgYZ - an array of real, x dimension is entirely local within 
! processors memory  while y and z are block-distributed among processors 
! in 2D grid
!
! The arrays may occupy the same memory space
! In this case their first elements should coincide. 
! Naturally, output overwrites input
!

      subroutine p3dfft_btran_c2r (XYZg,XgYZ)
!========================================================

      use fft_spec
      implicit none

      real(mytype),TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif

      integer x,y,z,i,k,nx,ny,nz
      integer(8) Nl

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif
      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

! FFT Tranform (C2C) in Z for all x and y 

      if(jproc .gt. 1) then

#ifdef STRIDE1
         call init_b_c(buf, 1,nz, &
              buf, 1, nz,nz,NB*iisize)
         call bcomm1_trans(XYZg,buf2,buf,timers(3),timers(10))
#else

         if(OW) then

            if(iisize*jjsize .gt. 0) then
               call init_b_c(XYZg, iisize*jjsize, 1, &
                    XYZg, iisize*jjsize, 1,nz,iisize*jjsize)

               timers(10) = timers(10) - MPI_Wtime()
               call exec_b_c2(XYZg, iisize*jjsize,1,  &
                    XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(10) = timers(10) + MPI_Wtime()
            endif
            call bcomm1(XYZg,buf,timers(3),timers(10))
   
         else
            if(iisize*jjsize .gt. 0) then
               call init_b_c(buf, iisize*jjsize, 1, &
                 buf, iisize*jjsize, 1,nz,iisize*jjsize)
            
               timers(10) = timers(10) - MPI_Wtime()
               call exec_b_c2(XYZg, iisize*jjsize,1, &
                 buf, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(10) = timers(10) + MPI_Wtime()
               call bcomm1(buf,buf,timers(3),timers(10))
            endif
         endif

#endif


      else
#ifdef STRIDE1
         call reorder_trans_b(XYZg,buf,buf2)
#else
         if(OW) then   
            Nl = iisize*jjsize*nz
            call init_b_c(XYZg, iisize*jjsize, 1, &
              XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
            call exec_b_c2(XYZg, iisize*jjsize, 1, &
              XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
            call ar_copy(XYZg,buf,Nl)
         else
            call init_b_c(XYZg, iisize*jjsize, 1, &
              buf, iisize*jjsize, 1,nz,iisize*jjsize)
            call exec_b_c2(XYZg, iisize*jjsize, 1, &
              buf, iisize*jjsize, 1,nz,iisize*jjsize)
         endif
#endif
      endif

! Exhange in columns if needed

!
! FFT Transform (C2C) in y dimension for all x, one z-plane at a time
!
      if(iisize * kjsize .gt. 0) then

         call init_b_c(buf,iisize,1,buf,iisize,1,ny,iisize)
         
         timers(9) = timers(9) - MPI_Wtime()
         do z=kjstart,kjend
               
            call btran_y_zplane(buf,z-kjstart,iisize,kjsize,iisize,1, &
                 buf,z-kjstart,iisize,kjsize,iisize,1,ny,iisize) 
            
         enddo
         timers(9) = timers(9) + MPI_Wtime()
      endif

      if(iproc .gt. 1) call bcomm2(buf,buf,timers(4),timers(10))

! Perform Complex-to-real FFT in x dimension for all y and z
      if(jisize * kjsize .gt. 0) then

         call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(8) = timers(8) - MPI_Wtime()
         call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(8) = timers(8) + MPI_Wtime()
      endif

      return
      end subroutine

      subroutine wrap_exec_b_c2(A,strideA,B,strideB, &
           N,m,L,k)

      complex(mytype) A(L,N),B(L,N)
      integer strideA,strideB,N,m,L,k
 

      call exec_b_c2(A(k,1),strideA,1,B(k,1),strideB,1,N,m)

      return
      end subroutine
