
! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with a 3D sine wave, then 
! performs forward transform, backward transform, and checks that 
! the results are correct, namely the same as in the start except 
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions. 
!
! Processor grid in this test is chosen to be square or close to square.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition, which 
! is the best in case it's feasible. 
!
! To build this program, use one of the provided makefiles, modifying 
! it as needed for your compilers, flags and library locations. 
! Be sure to link the program with both the P3DFFT library and the underlying 
! 1D FFT library such as ESSL or FFTW. 
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program ghost3d

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,n,m,x,y,z
      integer fstatus

      complex(mytype), dimension(:,:,:), allocatable, target::mem

      integer(8) Ntot
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer iproc,jproc

      integer nx,ny,nz
      common /grid/ nx,ny,nz

      integer gsize
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      common /p3dfftsets/ istart,iend,isize,fstart,fend,fsize,gsize

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      gsize = 2

      if (proc_id.eq.0) then 
         open (unit=3,file='stdin',status='old', access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif 
         ndim = 2

        read (3,*) nx, ny, nz, ndim,n
        write (*,*) "procs=",nproc," nx=",nx, " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
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
         dims(1) = 0
         dims(2) = 0
         call MPI_Dims_create(nproc,2,dims,ierr)
         if(dims(1) .gt. dims(2)) then
            dims(1) = dims(2)
            dims(2) = nproc / dims(1)
         endif
      endif

      iproc = dims(1)
      jproc = dims(2)

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif

      call p3dfft_setup (dims,nx,ny,nz,.false.)
      call p3dfft_init_ghosts(gsize)
      call get_dims(istart,iend,isize,1)
      call get_dims(fstart,fend,fsize,2)

      print *,'Allocating BEG (',isize,istart,iend
      allocate(mem(fsize(1)+2*gsize,fsize(2)+2*gsize,fsize(3)+2*gsize), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array mem3d'
      endif

      do i=1,n
      	call ghost_test(mem,mem)
      enddo

      call MPI_FINALIZE (ierr)

      stop

      end

!=========================================================

!     !========================================
!     !
!     !		ghost-cell test
!     !
!     !========================================
      subroutine ghost_test(mem_3d,mem_1d)
      	use p3dfft
      	use mpi
      	implicit none

!    !	common block
      	integer gsize
      	integer istart(3),iend(3),isize(3)
      	integer fstart(3),fend(3),fsize(3)
      	common /p3dfftsets/istart,iend,isize,fstart,fend,fsize,gsize

!     !	function args
      	real(mytype), dimension(istart(1):iend(1),istart(2):iend(2), &
                               istart(3):iend(3)), target	::	mem_3d
      	real(mytype), dimension(*), target	::	mem_1d

      	double precision  :: delta
      	double precision  :: maxdelta1,maxdelta2,maxdelta3
      	double precision  :: maxdelta01,maxdelta02,maxdelta03
      	integer            :: ierr, proc_id

!     	do ghost-cells check
      	call pencil_fill_3dsinus(mem_3d, ierr)
     	call update_rghosts(mem_3d)
     	call pencil_check_3dsinus_mneighbs(mem_1d, &
                  maxdelta1, maxdelta2, maxdelta3, ierr)

!     	print error
      	call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)
      	call MPI_reduce(maxdelta1,maxdelta01,1,MPI_DOUBLE_PRECISION,&
                 MPI_MAX,0,mpi_comm_world,ierr)
      	call MPI_reduce(maxdelta2,maxdelta02,1,MPI_DOUBLE_PRECISION,&
                 MPI_MAX,0,mpi_comm_world,ierr)
      	call MPI_reduce(maxdelta3,maxdelta03,1,MPI_DOUBLE_PRECISION,&
                 MPI_MAX,0,mpi_comm_world,ierr)
     	if(proc_id .eq. 0) then
      		if(maxdelta01 .gt. maxdelta02) maxdelta02 = maxdelta01
     		if(maxdelta02 .gt. maxdelta03) maxdelta03 = maxdelta02
     		write(*,*) '    Domain decomposition error   : ' &
     		            ,maxdelta03
     	endif

      return
      end subroutine ghost_test

!     !========================================
!     !
!     !		fill pencil with 3d-sinus-wave
!     !
!     !========================================
      subroutine pencil_fill_3dsinus(mem,ierr)
      use p3dfft
      implicit none

!    !	common block
      	integer nx,ny,nz
      	common /grid/ nx,ny,nz

      	integer gsize
      	integer istart(3),iend(3),isize(3)
      	integer fstart(3),fend(3),fsize(3)
      	common /p3dfftsets/ istart,iend, isize,fstart,fend,fsize,gsize

!     !	function args
      	real(mytype), dimension(istart(1):iend(1),istart(2):iend(2),&
                                istart(3):iend(3)), intent(out)::mem
      	integer, intent(out) :: ierr

!     !	other vars
      	integer :: i,j,k
      	double precision :: sinyz,pi2

      	pi2 = 8.d0 *datan(1.d0)
      	ierr = 0

!     !	fill array (pencil) with sinus-wave
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      sinyz = (sin((j-1)*pi2/ny)) *(sin((k-1)*pi2/nz))
      	      do i=istart(1),iend(1)
      	         mem(i,j,k)=(sin((i-1)*pi2/nx)) *sinyz
      	      enddo
      	   enddo
      	enddo

      end subroutine pencil_fill_3dsinus

!     !========================================
!     !
!     !		check pencil with 3d-sinus-wave
!     !
!     !========================================
      subroutine pencil_check_3dsinus_mneighbs(mem_1d, &
                                    maxdelta1,maxdelta2,maxdelta3,ierr)
      	use p3dfft
      	implicit none

!    !	common block
      	integer nx,ny,nz
      	common /grid/ nx,ny,nz

      	integer gsize
      	integer istart(3),iend(3),isize(3)
      	integer fstart(3),fend(3),fsize(3)
      	common /p3dfftsets/ istart,iend,isize,fstart,fend,fsize,gsize

!     !	function args
      	real(mytype), dimension(*) :: mem_1d
      	double precision, intent(out)::maxdelta1,maxdelta2,maxdelta3
      	integer, intent(out) :: ierr

!     !	other vars
      	integer :: i,j,k
      	double precision :: sin3_im, sin3_jm, sin3_km
      	double precision :: delta1, delta2, delta3
      	double precision :: pi2

      	pi2 = 8.d0 *datan(1.d0)
      	ierr = 0
      	maxdelta1 = 0.d0
      	maxdelta2 = 0.d0
      	maxdelta3 = 0.d0

!     !	check array (pencil) against sinus-wave
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_im = (sin((i-1-1)*pi2/nx)) *(sin((j  -1)*pi2/ny))&
                         *(sin((k  -1)*pi2/nz))
      	        sin3_jm = (sin((i  -1)*pi2/nx)) *(sin((j-1-1)*pi2/ny))&
                         *(sin((k  -1)*pi2/nz))
      	        sin3_km = (sin((i  -1)*pi2/nx)) *(sin((j  -1)*pi2/ny))&
                         *(sin((k-1-1)*pi2/nz))
      	        delta1 = abs( mem_1d(gr_ijk2i(i-1,j,k)) -sin3_im )
      	        delta2 = abs( mem_1d(gr_ijk2i(i,j-1,k)) -sin3_jm )
      	        delta3 = abs( mem_1d(gr_ijk2i(i,j,k-1)) -sin3_km )
      	        if( delta1 .gt. maxdelta1) then
      	        	maxdelta1 = delta1
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	        	maxdelta2 = delta2
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	        	maxdelta3 = delta3
      	        endif
      	      enddo
      	   enddo
      	enddo

      end subroutine pencil_check_3dsinus_mneighbs
