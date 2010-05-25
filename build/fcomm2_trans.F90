!========================================================
! Transpose an array from Y-pencils into Z-pencils
! and transform FFT in Z
!
! Assume Stride1 data structure
! Uses MPI_Alltoall(v)
!
      subroutine fcomm2_trans(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,ny_fft,kjsize)
!      complex(mytype) buf(nz_fft*iisize*jjsize)
      complex(mytype) dest(nz_fft,iisize,jjsize)

      real(8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,ix,x2,n
      integer*8 position,pos1

! Pack send buffers for exchanging y and z for all x at once 

      position = 1
 
      do i=0,jproc-1
         do z=1,kjsize
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo

#ifdef USE_EVEN
         position = position + (KfCntMax/(mytype*2) - jjsz(i)*iisize*kjsize)
#endif
      enddo
      
! Exchange y-z buffers in columns of processors

      t = t - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

      call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf,KfCntMax, mpi_byte,mpi_comm_col,ierr)

      t = MPI_Wtime() + t
      tc = tc - MPI_Wtime()
         
      do y=1,jjsize,NB
         y2 = min(y+NB-1,jjsize)

         do x=1,iisize,NB1
            x2 = min(x+NB1-1,iisize)
            
            n = iisize - (x2-x)
            do i=0,jproc-1
               pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
               do z=kjst(i),kjen(i)
                  position = pos1
                  do iy=y,y2
                     do ix=x,x2
! Here we are sure that dest and buf are different
                        dest(z,ix,iy) = buf(position)
                        position = position +1
                     enddo
                     position = position + n
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
         enddo
         call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
              nz_fft,NB*iisize)

      enddo
      tc = tc + MPI_Wtime()

#else
! Use MPI_Alltoallv

      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte, &
           buf,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = MPI_Wtime() + t

      tc = tc - MPI_Wtime()

      if(NB .eq. 1) then

         do y=1,jjsize

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)
            
               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do ix=x,x2
! Here we are sure that dest and buf are different
                     dest(z,ix,y) = buf(position)
                     position = position +1
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
            call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
                 nz_fft,iisize)
         enddo

      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)
               
               n = iisize - (x2-x)
               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do iy=y,y2
!                  position = ((z-1)*jjsize +iy-1)*iisize + x
                     do ix=x,x2
! Here we are sure that dest and buf are different
                        dest(z,ix,iy) = buf(position)
                        position = position +1
                     enddo
                     position = position + n
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
            call exec_f_c2(dest(1,1,y),1,nz_fft,dest(1,1,y),1,nz_fft, &
                 nz_fft,NB*iisize)
         enddo
         
      endif
      tc = tc + MPI_Wtime()
         
#endif
      return
      end subroutine
