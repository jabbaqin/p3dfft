!========================================================
! Transpose an array in Y=pencils into X-pencils
!
      subroutine bcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) dest(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(mytype) source(ny_fft,iisize,kjsize)
#else
      complex(mytype) source(iisize,ny_fft,kjsize)
#endif
      real(8) t,tc
      integer x,y,z,i,ierr,ix,iy,x2,y2,l
      integer OMP_GET_THREAD_NUM,ithr,st_y,en_y,sz_y
      integer(8) position,pos1,pos0,pos2

      
#ifdef OPENMP
      ithr = OMP_GET_THREAD_NUM()
#else
      ithr = 0
#endif

! !$OMP BARRIER
!$OMP MASTER 
      tc = tc - MPI_Wtime()
!$OMP END MASTER 

! c Pack and exchange x-z buffers in rows
      
      do i=0,iproc-1

      l = mod(jisz(i),num_thr)
      if(ithr .lt. num_thr-l) then
         sz_y = jisz(i)/num_thr
         st_y = jist(i) + sz_y*ithr
         en_y = st_y + sz_y-1
      else
         sz_y = jisz(i)/num_thr+1
         st_y = jist(i) + (jisz(i)/num_thr)*(num_thr-l) + sz_y*(ithr-num_thr+l)
         en_y = st_y + sz_y-1
      endif
      if(ithr .eq. num_thr-1) then
         en_y = jien(i)
      endif


#ifdef USE_EVEN
         pos0 = i*IfCntMax/(mytype*2)+1+ (st_y-jist(i))*iisize
#else
         pos0 = IfRcvStrt(i)/(mytype*2)+1+ (st_y-jist(i))*iisize
#endif

         do z=1,kjsize
            
#ifdef STRIDE1
            pos1 = pos0
            do y=st_y,en_y,nby1
               y2 = min(y+nby1-1,en_y)
               do x=1,iisize,nbx
                  x2 = min(x+nbx-1,iisize)
                  pos2 = pos1 + x-1 
                  do iy = y,y2
                     position = pos2
                     do ix=x,x2
                        buf1(position) = source(iy,ix,z)
                        position = position + 1
                     enddo
                     pos2 = pos2 + iisize
                  enddo
               enddo
               pos1 = pos1 + iisize*nby1
            enddo

#else
            position = pos0 
            do y=st_y,en_y
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position + 1
               enddo
            enddo
#endif

            pos0 = pos0 + iisize*jisz(i)
         enddo
      enddo

!$OMP BARRIER
!$OMP MASTER
      tc = tc + MPI_Wtime() 
      t = t - MPI_Wtime() 

#ifdef USE_EVEN
      call mpi_alltoall (buf1,IfCntMax, mpi_byte, buf2,IfCntMax,mpi_byte, mpi_comm_row,ierr)
#else
      call mpi_alltoallv (buf1,KrSndCnts, KrSndStrt, mpi_byte, buf2,KrRcvCnts,KrRcvStrt,mpi_byte, mpi_comm_row,ierr)
#endif

      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()
!$OMP END MASTER
!$OMP BARRIER

! Unpack receive buffers into dest

      l = mod(jisize,num_thr)
      if(ithr .lt. num_thr-l) then
         sz_y = jisize/num_thr
         st_y = 1 + sz_y*ithr
         en_y = st_y + sz_y-1
      else
         sz_y = jisize/num_thr+1
         st_y = 1 + (jisize/num_thr)*(num_thr-l) + sz_y*(ithr-num_thr+l)
         en_y = st_y + sz_y-1
      endif
      if(ithr .eq. num_thr-1) then
         en_y = jisize
      endif

      do i=0,iproc-1
#ifdef USE_EVEN
         pos0 = i*IfCntMax/(mytype*2) + 1 + (st_y-1)*iisz(i) 
#else
         pos0 = IfSndStrt(i)/(mytype*2) + 1 + (st_y-1)*iisz(i) 
#endif
         do z=1,kjsize
            position = pos0 
            do y=st_y,en_y
               do x=iist(i),iien(i)
                  dest(x,y,z) = buf2(position)
                  position = position +1
               enddo
            enddo
            pos0 = pos0 + iisz(i)*jisize
         enddo
      enddo

!$OMP BARRIER
!$OMP MASTER 
      tc = tc + MPI_Wtime() 
!$OMP END MASTER 
      
      return
      end subroutine
