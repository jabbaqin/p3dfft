!========================================================
! Transpose array in X into Y-pencils
! Uses MPI_Alltoall(v)
! 
      subroutine fcomm1(source,dest,t)
!========================================================

      implicit none

      complex(mytype) source(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(mytype) dest(ny_fft,iisize,kjsize)
#else
      complex(mytype) dest(iisize,ny_fft,kjsize)
#endif

      real(8) t
      integer x,y,i,ierr,z,xs,j,n,ix,iy
      integer(8) position,pos1

! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf

      position = 1

      do i=0,iproc-1
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  buf1(position) = source(x,y,z)
                  position = position +1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1)*IfCntMax/(mytype*2)+1
#endif         
      enddo
      
#ifdef USE_EVEN

! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors) 

      t = t - MPI_Wtime()
      call mpi_alltoall(buf1,IfCntMax, mpi_byte, &
           buf2,IfCntMax, mpi_byte,mpi_comm_row,ierr)
      t = t + MPI_Wtime()

#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
      t = t - MPI_Wtime() 
      call mpi_alltoallv(buf1,IfSndCnts, IfSndStrt,mpi_byte, &
           buf2,IfRcvCnts, IfRcvStrt,mpi_byte,mpi_comm_row,ierr)
      t = MPI_Wtime() + t
#endif

! Unpack the data

#ifdef STRIDE1
      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            pos1 = position
            do y=jist(i),jien(i),nb1
               do x=1,iisize,nb1
                  do iy = y,min(y+nb1-1,jien(i))
                     position = pos1 + x-1 +(iy-jist(i))*iisize
                     do ix=x,min(x+nb1-1,iisize)
                        dest(iy,ix,z) = buf2(position)
                        position = position + 1
                     enddo
                  enddo
               enddo
            enddo
            position = pos1 + jisz(i) * iisize
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo

#else
      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            do y=jist(i),jien(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position + 1
               enddo
            enddo
         enddo
#ifdef USE_EVEN
         position = (i+1) * IfCntMax/(mytype*2) + 1
#endif
      enddo
#endif

      return
      end subroutine
