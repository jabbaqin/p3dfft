! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
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

      integer i,n,nx,ny,nz,nv,j
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real(mytype), dimension(:),  allocatable :: B
      real(mytype) pi,twopi,sinyz,diff,cdiff,ccdiff,ans

      integer(i8) Ntot
      real(mytype) factor
      real(mytype),dimension(:,:),allocatable:: sinx,siny,sinz
      real(r8) rtime1,rtime2,Nglob,prec
      real(r8) gt(12,3),gtcomm(3),tc
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer, dimension(3) :: istart,iend,isize
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

        read (3,*) nx, ny, nz, ndim,nv,n
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"ndim=",ndim,"nv=",nv," repeat=", n
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
      call MPI_Bcast(nv,1, MPI_INTEGER,0,mpi_comm_world,ierr)
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

!      allocate (B(istart(1):istart(1)+memsize(1)-1,istart(2):istart(2)+memsize(2)-1,istart(3):istart(3)+memsize(3)-1,nv), stat=ierr)

      allocate (B(memsize(1)*memsize(2)*memsize(3)*nv), stat=ierr)

      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array B'
      endif


! Initialize the array to be transformed
!
	call init_ar(B,nv)
!
! transform from physical space to wavenumber space
! (XgYiZj to XiYjZg)
! then transform back to physical space
! (XiYjZg to XgYiZj)
!
! Time the calls with MPI_Wtime

         call ftran_r2c_many (B,isize(1)*isize(2)*isize(3),B, &
             fsize(1)*fsize(2)*fsize(3),nv,'fft')

         call ftran_r2c_many (B,isize(1)*isize(2)*isize(3),B, &
             fsize(1)*fsize(2)*fsize(3),nv,'fft')

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

         call ftran_r2c_many (B,isize(1)*isize(2)*isize(3),B, &
             fsize(1)*fsize(2)*fsize(3),nv,'fft')
         
         rtime1 = rtime1 + MPI_wtime()

           call print_all_real(B,nv,proc_id,Nglob)
         
! Normalize
           call mult_array(B, nv,Ntot*2,factor)

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()

! Backward transform: note - we pass the same array both as 
! original (complex) and transformed (real)

         call btran_c2r_many (B,fsize(1)*fsize(2)*fsize(3),B, &
             isize(1)*isize(2)*isize(3),nv,'tff')
         rtime1 = rtime1 + MPI_wtime()
         
      end do

! Clean the FFT work space
      call p3dfft_clean

! Check results

     call check_res(B,nv,Nglob)

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

      subroutine check_res(B,nv,Nglob)

      real(mytype) B(isize(1),isize(2),isize(3),nv)
      integer nv,j,x,y,z
!      integer, dimension(3) :: isize
      real(mytype) cdiff,sinyz,ans,ccdiff
      real(r8) prec,Nglob

      do j=1,nv

         cdiff=0.0d0
         do 20 z=1,isize(3)
            do 20 y=1,isize(2)
               sinyz=siny(y,j)*sinz(z,j)
               do 20 x=1,isize(1)
                  ans=sinx(x,j)*sinyz
                  if(cdiff .lt. abs(B(x,y,z,j)-ans)) then
                     cdiff = abs(B(x,y,z,j)-ans)
                  endif
 20      continue
         call MPI_Reduce(cdiff,ccdiff,1,mpireal,MPI_MAX,0, &
                   MPI_COMM_WORLD,ierr)

         if (proc_id.eq.0) write (6,*) 'Var. ',j,': max diff =',ccdiff

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
      enddo

      return
      end subroutine


      subroutine init_ar(B,nv)
      
      real(mytype) B(isize(1),isize(2),isize(3),nv),sinyz
      integer nv,j,x,y,z
!      integer, dimension(3) :: isize

      allocate (sinx(nx,nv))
      allocate (siny(ny,nv))
      allocate (sinz(nz,nv))

      do j=1,nv

         do z=istart(3),iend(3)
            sinz(z-istart(3)+1,j)=sin(j*(z-1)*twopi/nz)
         enddo
         do y=istart(2),iend(2)
            siny(y-istart(2)+1,j)=sin(j*(y-1)*twopi/ny)
         enddo
         do x=istart(1),iend(1)
            sinx(x-istart(1)+1,j)=sin(j*(x-1)*twopi/nx)
         enddo
         
         do z=1,isize(3)
            do y=1,isize(2)
               sinyz=siny(y,j)*sinz(z,j)
               do x=1,isize(1)
                  B(x,y,z,j)=sinx(x,j)*sinyz 
               enddo
            enddo
         enddo
      enddo

      return
      end subroutine

!=========================================================

      subroutine mult_array(X,nv,nar,f)

      use p3dfft

      integer(i8) nar,i
      integer nv,j
      real(mytype) X(nar,nv)
      real(mytype) f

      do j=1,nv
      do i=1,nar
         X(i,j) = X(i,j) * f
      enddo
      enddo

      return
      end subroutine

!=========================================================
! Translate one-dimensional index into three dimensions,
!    print out significantly non-zero values
!
      subroutine print_all_real(Ar,nv,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id,j,nv
      real(r8) Nglob
      real(mytype), target :: Ar(Fsize(1)*2,Fsize(2),Fsize(3),nv)
!      integer Fstart(3),Fend(3),Fsize(3)

!      call p3dfft_get_dims(Fstart,Fend,Fsize,2)

         do j=1,nv
           if(proc_id .eq. 0) then
              print *,'Result of forward transform, var. ',j
           endif

!      do i=1,Nar*2,2
!            z = (i-1)/(Fsize(1)*Fsize(2)*2)
!            y = ((i-1)/2 - z * Fsize(1)*Fsize(2))/(Fsize(1))
!            x = (i-1)/2-z*Fsize(1)*Fsize(2) - y*Fsize(1)

	     do z=1,Fsize(3)
	     do y=1,Fsize(2)
	     do x=1,Fsize(1)*2,2
         if(abs(Ar(x,y,z,j)) + abs(Ar(x+1,y,z,j)) .gt. Nglob *1.25e-6) then
            print *,'(',x+Fstart(1)-1,y+Fstart(2)-1,z+Fstart(3)-1,') ',Ar(x,y,z,j),Ar(x+1,y,z,j)
         endif
	 enddo
	 enddo
      enddo
      enddo

      return
      end subroutine


      end
