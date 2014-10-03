! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
!    Copyright (C) 2010-2011 Jens Henrik Goebbert
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

!========================================================
! Transpose Y and Z pencils
! Assume Stride1 data structure

      subroutine fcomm2_trans(source,dest,op,t,tc)
!========================================================

      use fft_spec
      implicit none

! Assume stride1
      integer nz,dim_out
      complex(mytype) source(ny_fft,iisize,kjsize)
      complex(mytype) dest(nzc,jjsize,iisize)
!      complex(mytype) buf3(nz_fft,jjsize)
      complex(mytype), allocatable :: buf1(:),buf2(:)

      real(r8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,ix,x2,n,sz,l,dny,dnz
      integer(i8) position,pos1,pos0,pos2
      integer threadid,omp_get_thread_num,ithr
      character(len=3) op
      integer sndcnts(0:jproc-1)
      integer rcvcnts(0:jproc-1)
      integer sndstrt(0:jproc-1)
      integer rcvstrt(0:jproc-1)


! Pack send buffers for exchanging y and z for all x at once 


!      tc = tc - MPI_Wtime()

!      if(taskid .eq. 0) then
!         print *,'Entering fcomm2_trans; source='
!      endif
!      do j=1,nv
!      print *,taskid,':, j=',j
!      call print_buf(source(1,1,1,j),ny_fft,iisize,kjsize)
!      enddo


     allocate(buf1(buf_size))
     allocate(buf2(buf_size))

      threadid = OMP_GET_THREAD_NUM()
#ifdef DEBUG
      print *,taskid,threadid,": Entered fcomm2_trans"
!$OMP FLUSH
#endif

     call pack_fcomm2_trans(buf1,source)

!      tc = tc + MPI_Wtime()

#ifdef DEBUG
      print *,taskid,threadid,": Passed barrier in fcomm2_trans"
!$OMP FLUSH
#endif
!      t = t - MPI_Wtime()


!  can't use !$OMP ORDERED for the second time in the loop, so write our own

   do ithr=0,nthreads-1
!$OMP BARRIER
      if(ithr .eq. threadid) then   

#ifdef USE_EVEN
      call mpi_alltoall(buf1,KfCntMax, mpi_byte, buf2,KfCntMax, mpi_byte,mpi_comm_col,ierr)
#else      
! Exchange y-z buffers in columns of processors

#ifdef DEBUG
      print *,taskid,threadid,": Calling alltoall in fcomm2_trans"
!$OMP FLUSH
#endif
      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte,buf2,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif

	endif
    enddo	

!  !$OMP END ORDERED

#ifdef DEBUG
      print *,taskid,threadid,": Passed alltoall in fcomm2_trans"
!$OMP FLUSH
#endif

!      t = MPI_Wtime() + t

	deallocate(buf1)

     if(jjsize .gt. 0) then

!      tc = - MPI_Wtime() + tc
	
      call unpack_fcomm2_trans(dest,buf2,op)

!      tc = tc + MPI_Wtime()

     endif


!      if(taskid .eq. 0) then
!         print *,'Exiting fcomm2_trans; dest='
!      endif
!      do j=1,nv
!      print *,taskid,':, j=',j
!      call print_buf(dest(1,1,1,j),nz_fft,jjsize,iisize)
!      enddo
	

	deallocate(buf2)
!	enddo	

#ifdef DEBUG
      print *,taskid,threadid,": Exiting fcomm2_trans"
!$OMP FLUSH
#endif

      return
      end subroutine

!------------------------------------------------------------
      subroutine pack_fcomm2_trans(sendbuf,source)

      use fft_spec
      implicit none

      complex(mytype) source(ny_fft,iisize,kjsize)
#ifdef USE_EVEN
      complex(mytype) sendbuf(KfCntMax)
#else
      complex(mytype) sendbuf(ny_fft*iisize*kjsize)
#endif
      integer x,z,y,i,y2,z2,iy,iz,ix,x2,dny
      integer(i8) position,pos1,pos0,pos2
      integer threadid,OMP_GET_THREAD_NUM


#ifdef DEBUG
      threadid = OMP_GET_THREAD_NUM()
      print *,taskid,threadid,": Entered pack_fcomm2_trans"
!$OMP FLUSH
#endif
      dny = ny_fft-nyc

      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i * KfCntMax/(mytype*2)  + 1 
#else
         pos0 = KfSndStrt(i) /(mytype*2)+ 1 
#endif


! Pack the sendbuf, omitting the center ny-nyc elements in Y dimension

! If clearly in the first half of ny
         if(jjen(i) .le. nyhc) then
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do x=1,iisize
                  do y=jjst(i),jjen(i)
                     sendbuf(position) = source(y,x,z)
                     position = position+1
                  enddo
               enddo	
            enddo

! If clearly in the second half of ny
         else if (jjst(i) .ge. nyhc+1) then
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do x=1,iisize
                  do y=jjst(i)+dny,jjen(i)+dny
                     sendbuf(position) = source(y,x,z)
                     position = position+1
                  enddo
               enddo	
            enddo



! If spanning the first and second half of ny (e.g. iproc is odd)
         else
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do x=1,iisize
                  do y=jjst(i),nyhc
                     sendbuf(position) = source(y,x,z)
                     position = position+1
                  enddo
                  do y=ny_fft-nyhc+1,jjen(i)+dny
                     sendbuf(position) = source(y,x,z)
                     position = position+1
                  enddo
               enddo	
            enddo
         endif

      enddo

#ifdef DEBUG
      print *,taskid,threadid,": Exiting pack_fcomm2_trans"
!$OMP FLUSH
#endif
      return
      end subroutine

!------------------------------------------------------------
      subroutine unpack_fcomm2_trans(dest,recvbuf,op)

      use fft_spec
      implicit none

      complex(mytype) dest(nzc,jjsize,iisize)
      complex(mytype) buf3(nz_fft,jjsize)
#ifdef USE_EVEN
      complex(mytype) recvbuf(KfCntMax)
#else
      complex(mytype) recvbuf(nzc*jjsize*iisize)
#endif
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,ix,x2,n,sz,l,dny,dnz,nz,dim_out
      integer(i8) position,pos1,pos0,pos2
      character(len=3) op

      nz = nz_fft
      if(nz .ne. nzc) then
	dnz = nz - nzc


      if(op(3:3) == '0' .or. op(3:3) == 'n') then

      do x=1,iisize
	pos0 = (x-1)*jjsize

         do i=0,jproc-1

#ifdef USE_EVEN
 	      pos1 = pos0 + i  * KfCntMax / (mytype*2) 
#else
 	      pos1 = pos0 + KfRcvStrt(i) / (mytype*2) 
#endif

           if(kjen(i) .lt. nzhc .or. kjst(i) .gt. nzhc+1) then
! Just copy the data
              do z=kjst(i),kjen(i),NBz
                 z2 = min(z+NBz-1,kjen(i))

                 do y=1,jjsize,NBy2
                    y2 = min(y+NBy2-1,jjsize)
                    pos2 = pos1 + y 

                    do iz=z,z2
                       position = pos2
                       do iy=y,y2
! Here we are sure that dest and buf are different
                         dest(iz,iy,x) = recvbuf(position)
                         position = position +1
                       enddo
                       pos2 = pos2 + iisize * jjsize
                    enddo
                 enddo
                 pos1 = pos1 + iisize*jjsize*NBz
              enddo

	    else 
! Copy some data, then insert zeros to restore full dimension, then again 
! copy some data if needed	
               do z=kjst(i),nzhc,NBz
       	          z2 = min(z+NBz-1,nzhc)

                 do y=1,jjsize,NBy2
                    y2 = min(y+NBy2-1,jjsize)
                    pos2 = pos1 + y 

                    do iz=z,z2
                       position = pos2
                       do iy=y,y2
! Here we are sure that dest and buf are different
                         dest(iz,iy,x) = recvbuf(position)
                         position = position +1
                       enddo
                       pos2 = pos2 + iisize * jjsize
                    enddo
                 enddo
                 pos1 = pos1 + iisize*jjsize*NBz
              enddo

	      pos0 = pos0 + iisize*jjsize*dnz
#ifdef USE_EVEN
 	      pos1 = pos0 + i * KfCntMax / (mytype*2) 
#else
 	      pos1 = pos0 + KfRcvStrt(i) / (mytype*2) 
#endif

	      do z=nzhc+1,kjen(i),NBz
                 z2 = min(z+NBz-1,kjen(i))

                 do y=1,jjsize,NBy2
                    y2 = min(y+NBy2-1,jjsize)
                    pos2 = pos1 + y 

                    do iz=z,z2
                       position = pos2
                       do iy=y,y2
! Here we are sure that dest and buf are different
                         dest(iz,iy,x) = recvbuf(position)
                         position = position +1
                       enddo
                       pos2 = pos2 + iisize * jjsize
                    enddo
                 enddo
                 pos1 = pos1 + iisize*jjsize*NBz
              enddo
	   endif

	enddo

       enddo	

     else
 
       do x=1,iisize

         do i=0,jproc-1

#ifdef USE_EVEN
            pos0 = i  *KfCntMax / (mytype*2) + (x-1)*jjsize 
#else         
	    pos0 = KfRcvStrt(i) / (mytype*2) + (x-1)*jjsize
#endif 

            do z=kjst(i),kjen(i),NBz
               z2 = min(z+NBz-1,kjen(i))
            
               do y=1,jjsize,NBy2
                  y2 = min(y+NBy2-1,jjsize)

                  pos1 = pos0 +y

                  do iz=z,z2
                     position = pos1 
                     do iy=y,y2
! Here we are sure that dest and buf are different
                        buf3(iz,iy) = recvbuf(position)
                        position = position + 1
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
               pos0 = pos0 + jjsize*iisize*NBz
           enddo 
       enddo	

	if(op(3:3) == 't' .or. op(3:3) == 'f') then
             call exec_f_c2_same(buf3, 1,nz_fft, &
			  buf3, 1,nz_fft,nz_fft,jjsize)
	else if(op(3:3) == 'c') then	
             call exec_ctrans_r2_complex_same(buf3, 2,2*nz_fft, &
			  buf3, 2,2*nz_fft,nz_fft,jjsize)
	else if(op(3:3) == 's') then	
             call exec_strans_r2_complex_same(buf3, 2,2*nz_fft, &
		          buf3, 2,2*nz_fft,nz_fft,jjsize)
 	else
	   print *,taskid,'Unknown transform type: ',op(3:3)
	   call MPI_abort(MPI_COMM_WORLD,ierr)
	endif

	do y=1,jjsize	
	   do z=1,nzhc
	      dest(z,y,x) = buf3(z,y)
	   enddo
	   do z=nzhc+1,nzc
	      dest(z,y,x) = buf3(z+dnz,y)
	   enddo
	enddo

      enddo
      endif
      else  ! if nz = nzc

      if(op(3:3) == '0' .or. op(3:3) == 'n') then

      do x=1,iisize

         do i=0,jproc-1

#ifdef USE_EVEN
            pos0 = i *KfCntMax / (mytype*2) + (x-1)*jjsize 
#else         
	   pos0 = KfRcvStrt(i) / (mytype*2) + (x-1)*jjsize
#endif 

            do z=kjst(i),kjen(i),NBz
               z2 = min(z+NBz-1,kjen(i))
            
               do y=1,jjsize,NBy2
                  y2 = min(y+NBy2-1,jjsize)

                  pos1 = pos0 +y

                  do iz=z,z2
                     position = pos1 
                     do iy=y,y2
! Here we are sure that dest and buf are different
                        dest(iz,iy,x) = recvbuf(position)
                        position = position + 1
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
               pos0 = pos0 + jjsize*iisize*NBz
           enddo 
        enddo
      enddo

      else
      
      do x=1,iisize

         do i=0,jproc-1

#ifdef USE_EVEN
            pos0 = i *KfCntMax / (mytype*2) + (x-1)*jjsize 
#else         
	   pos0 = KfRcvStrt(i) / (mytype*2) + (x-1)*jjsize
#endif 

            do z=kjst(i),kjen(i),NBz
               z2 = min(z+NBz-1,kjen(i))
            
               do y=1,jjsize,NBy2
                  y2 = min(y+NBy2-1,jjsize)

                  pos1 = pos0 +y

                  do iz=z,z2
                     position = pos1 
                     do iy=y,y2
! Here we are sure that dest and buf are different
                        buf3(iz,iy) = recvbuf(position)
                        position = position + 1
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
               pos0 = pos0 + jjsize*iisize*NBz
            enddo
         enddo

	if(op(3:3) == 't' .or. op(3:3) == 'f') then
             call exec_f_c2_dif(buf3, 1,nz_fft, &
			  dest(1,1,x), 1,nz_fft,nz_fft,jjsize)
	else if(op(3:3) == 'c') then	
             call exec_ctrans_r2_complex_dif(buf3, 2,2*nz_fft, &
			  dest(1,1,x), 2,2*nz_fft,nz_fft,jjsize)
	else if(op(3:3) == 's') then	
             call exec_strans_r2_complex_dif(buf3, 2,2*nz_fft, &
		          dest(1,1,x), 2,2*nz_fft,nz_fft,jjsize)
 	else
	   print *,taskid,'Unknown transform type: ',op(3:3)
	   call MPI_abort(MPI_COMM_WORLD,ierr)
	endif

      enddo

      endif
      endif

      return
      end subroutine

