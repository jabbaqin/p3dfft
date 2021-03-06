! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2014 Dmitry Pekurovsky
!    Copyright (C) 2006-2014 University of California
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!----------------------------------------------------------------------------

! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with a 3D sine wave, then 
! performs 3D FFT forward transform IN PLACE, then backward transform 
! In PLACE, and checks that 
! the results are correct, namely the same as in the start except 
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions. 
!
! The program expects 'stdin' file in the working directory, with 
! a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here Nx,Ny,Nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice 
! of processor geometry in case of 2D decomposition. It should contain 
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program fft3d

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real(mytype), dimension(:,:,:),  allocatable :: B
      real(mytype) twopi

      integer(i8) Ntot
      real(mytype) factor
      real(mytype),dimension(:),allocatable:: sinx,siny,sinz
      real(r8) rtime1,rtime2,Nglob
      real(r8) gt(12,3),gtcomm(3),tc
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      integer iproc,jproc,add
      logical iex
      integer memsize(3),tmp1(3),tmp2(3)

! Initialize MPI

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      twopi=atan(1.0d0)*8.0d0

      timers = 0.0
      gt=0.0
      gtcomm=0.0

      if (proc_id.eq.0) then 
         open (unit=3,file='stdin',status='old', &
               access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif 
         ndim = 2

        read (3,*) nx, ny, nz, ndim,n
	print *,'P3DFFT test, 3D wave input, in-place transform'
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
        if(mytype .eq. 4) then
           print *,'Single precision version'
        else if(mytype .eq. 8) then
           print *,'Double precision version'
        endif
       endif

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ndim,1, MPI_INTEGER,0,mpi_comm_world,ierr)

!    nproc is devided into a iproc x jproc stencle
!

      if(ndim .eq. 1) then
         dims(1) = 1
         dims(2) = nproc
      else if(ndim .eq. 2) then
         inquire(file='dims',exist=iex)
         if (iex) then
            if (proc_id.eq.0) print *, 'Reading proc. grid from file dims'
            open (999,file='dims')
            read (999,*) dims(1), dims(2)
            close (999)
            if(dims(1) * dims(2) .ne. nproc) then
               dims(2) = nproc / dims(1)
            endif
         else
            if (proc_id.eq.0) print *, 'Creating proc. grid with mpi_dims_create'
            dims(1) = 0
            dims(2) = 0
            call MPI_Dims_create(nproc,2,dims,ierr)
            if(dims(1) .gt. dims(2)) then
               dims(1) = dims(2)
               dims(2) = nproc / dims(1)
            endif
         endif
      endif
       
      iproc = dims(1)
      jproc = dims(2)

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif

! Set up work structures for P3DFFT
      call p3dfft_setup (dims,nx,ny,nz,MPI_COMM_WORLD,nx,ny,nz,.true.)

! Get dimensions for the original array of real numbers, X-pencils
      call p3dfft_get_dims(istart,iend,isize,1)

! Get dimensions for the R2C-forward-transformed array of complex numbers
!   Z-pencils (depending on how the library was compiled, the first 
!   dimension could be either X or Z)
! 
      call p3dfft_get_dims(fstart,fend,fsize,2)


! Since we are allocating the same array for input and output, 
! we need to make sure it has enough space. This is achieved
! by calling get_dims with option 3. 

      call p3dfft_get_dims(tmp1,tmp2,memsize,3)

!      print *,proc_id,': memsize=',memsize

      allocate (B(istart(1):istart(1)+memsize(1)-1,istart(2):istart(2)+memsize(2)-1,istart(3):istart(3)+memsize(3)-1), stat=ierr)

      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array B'
      endif


! Initialize the array to be transformed
!
	call init_ar_sine(B)
!
! transform from physical space to wavenumber space
! (XgYiZj to XiYjZg)
! then transform back to physical space
! (XiYjZg to XgYiZj)
!
! Time the calls with MPI_Wtime

      Ntot = fsize(1)*fsize(2)
      Ntot = Ntot * fsize(3)
      Nglob = nx * ny
      Nglob = Nglob * nz
      factor = 1.0d0/Nglob

      rtime1 = 0.0               
      do  m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif
         
! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()

! Forward transform: note - we pass the same array both as 
! original (real) and transformed (complex)

         call ftran_r2c (B,B,'fft')
         
         rtime1 = rtime1 + MPI_wtime()

         if(proc_id .eq. 0) then
            print *,'Result of forward transform:'
         endif
         call print_all_real(B,Ntot,proc_id,Nglob)
         
! Normalize
         call mult_array(B, Ntot*2,factor)
         
! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()

! Backward transform: note - we pass the same array both as 
! original (complex) and transformed (real)

         call btran_c2r (B,B,'tff')       
         rtime1 = rtime1 + MPI_wtime()
         
      end do

! Clean the FFT work space
      call p3dfft_clean

! Check results
      call check_res(B)

! Gather timing statistics

      timers = timers / dble(n)

      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop', &
         proc_id,rtime2/dble(n)

      call MPI_Reduce(timers,gt(1,1),12,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,2),12,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,3),12,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      tc = (timers(1)+timers(2)+timers(3)+timers(4))
      call MPI_Reduce(tc,gtcomm(1),1,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(2),1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(3),1,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      gt(1:12,1) = gt(1:12,1) / dble(nproc)
      gtcomm(1) = gtcomm(1) / dble(nproc)

      if(proc_id .eq. 0) then
         do i=1,12
            print *,'timer',i,' (avg/max/min): ',gt(i,:)
         enddo
         print *,'Total comm (avg/max/min): ',gtcomm
      endif


      call MPI_FINALIZE (ierr)

      contains 

!=========================================================
	subroutine check_res(C)
!=========================================================

	real(mytype) C(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
	real(mytype) cdiff,ccdiff,sinyz,ans,prec
	integer x,y,z


      cdiff=0.0d0
      do 20 z=istart(3),iend(3)
         do 20 y=istart(2),iend(2)
            sinyz=siny(y)*sinz(z)
            do 20 x=istart(1),iend(1)
            ans=sinx(x)*sinyz
            if(cdiff .lt. abs(C(x,y,z)-ans)) then
               cdiff = abs(C(x,y,z)-ans)
!               print *,'x,y,z,cdiff=',x,y,z,cdiff
            endif
 20   continue
      call MPI_Reduce(cdiff,ccdiff,1,mpireal,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if(proc_id .eq. 0) then
         if(mytype .eq. 8) then
            prec = 1e-14
         else
            prec = 1e-5
         endif
         if(ccdiff .gt. prec * Nglob*0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
         write (6,*) 'max diff =',ccdiff
      endif

      return
      end subroutine

!=========================================================
	subroutine init_ar_sine(A)
!=========================================================

	real(mytype) A(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
	integer x,y,z
	real(mytype) sinyz

      allocate (sinx(nx))
      allocate (siny(ny))
      allocate (sinz(nz))

      do z=istart(3),iend(3)
         sinz(z)=sin((z-1)*twopi/nz)
      enddo
      do y=istart(2),iend(2)
         siny(y)=sin((y-1)*twopi/ny)
      enddo
      do x=istart(1),iend(1)
         sinx(x)=sin((x-1)*twopi/nx)
      enddo

! Initialize with 3D sine wave

      do z=istart(3),iend(3)
         do y=istart(2),iend(2)
            sinyz=siny(y)*sinz(z)
            do x=istart(1),iend(1)
               A(x,y,z)=sinx(x)*sinyz 
            enddo
         enddo
      enddo

      return
      end subroutine

!=========================================================

      subroutine mult_array(X,nar,f)

      use p3dfft

      integer(i8) nar,i
      real(mytype) X(nar)
      real(mytype) f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
      end subroutine

!=========================================================
! Translate one-dimensional index into three dimensions,
!    print out significantly non-zero values
!
      subroutine print_all_real(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(i8) i,Nar
      real(r8) Nglob
      real(mytype), target :: Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)

      call p3dfft_get_dims(Fstart,Fend,Fsize,2)
      do i=1,Nar*2,2
         if(abs(Ar(1,1,i)) + abs(Ar(1,1,i+1)) .gt. Nglob *1.25e-6) then
            z = (i-1)/(Fsize(1)*Fsize(2)*2)
            y = ((i-1)/2 - z * Fsize(1)*Fsize(2))/(Fsize(1))
            x = (i-1)/2-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i),Ar(1,1,i+1)
         endif
      enddo

      return
      end subroutine


      end
