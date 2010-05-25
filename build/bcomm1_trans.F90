!========================================================
! Transpose an array in Z=pencils  into Y-pencils
! and perform FFT transform in Z dimension
!
! Assumes stride1 data structure

      subroutine bcomm1_trans (source,buf3,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(nz_fft,iisize,jjsize)
      complex(mytype) buf3(nz_fft,iisize,NB)
      complex(mytype) dest(iisize,ny_fft,kjsize)

      real(8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n
      integer*8 position,pos1
      
!     Pack the data for sending


#ifdef USE_EVEN
! Use MPI_Alltoall

      tc = tc - MPI_Wtime()

      if(OW) then

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

            call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
              nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)
 
               do i=0,jproc-1
                  pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
                  do z=kjst(i),kjen(i)
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = source(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + iisize - (x2-x)
                     enddo
                     pos1 = pos1 + iisize * jjsize 
                  enddo
               enddo
            enddo
         enddo

         tc = tc + MPI_Wtime()
         t = t - MPI_Wtime() 
         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
         t = t + MPI_Wtime() 

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
! Here we are sure that buf is not the same as dest
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo
         
      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)
            
            call exec_b_c2(source(1,1,y),1,nz_fft,buf3(1,1,y),1,nz_fft, &
                 nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)
 
               do i=0,jproc-1
                  pos1 = KfCntMax/(mytype*2) * i +(y-1)*iisize+x
                  do z=kjst(i),kjen(i)
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = buf3(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + iisize - (x2-x)
                     enddo
                     pos1 = pos1 + iisize * jjsize 
                  enddo
               enddo
            enddo

         enddo
         
         tc = tc + MPI_Wtime()
         t = t - MPI_Wtime() 
         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
         
         t = t + MPI_Wtime() 

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo

      endif
         
#else
! Use MPI_Alltoallv

      tc = tc - MPI_Wtime()

      if(OW) then

         if(NB .eq.1) then

            do y=1,jjsize

               call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
                  nz_fft,iisize)

               do x=1,iisize,NB1
                  x2 = min(x+NB1-1,iisize)
 
                  pos1 = (y-1)*iisize + x
                  do z=1,nz_fft
                     position = pos1
                     do ix=x,x2
                        buf1(position) = source(z,ix,y)
                        position = position+1
                     enddo
                     pos1 = pos1 + iisize*jjsize
                  enddo
               enddo
            enddo
            
         else

            do y=1,jjsize,NB
               y2 = min(y+NB-1,jjsize)
               
               call exec_b_c2(source(1,1,y),1,nz_fft,source(1,1,y),1,nz_fft, &
                    nz_fft,NB*iisize)
               
               do x=1,iisize,NB1
                  x2 = min(x+NB1-1,iisize)
                  
                  n = iisize - (x2-x)
                  pos1 = (y-1)*iisize + x
                  do z=1,nz_fft
                     position = pos1
                     do iy = y,y2
                        do ix=x,x2
                           buf1(position) = source(z,ix,iy)
                           position = position+1
                        enddo
                        position = position + n
                     enddo
                     pos1 = pos1 + iisize*jjsize
                  enddo
               enddo
            enddo
            
         endif
         
         tc = tc + MPI_Wtime()

         t = t - MPI_Wtime() 
         call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, &
              buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
         t = t + MPI_Wtime() 

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
! Here we are sure that buf is not the same as dest
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo

      else

         do y=1,jjsize,NB
            y2 = min(y+NB-1,jjsize)

            call exec_b_c2(source(1,1,y),1,nz_fft,buf3(1,1,y),1,nz_fft, &
              nz_fft,NB*iisize)

            do x=1,iisize,NB1
               x2 = min(x+NB1-1,iisize)
 
               n = iisize - (x2-x)
               pos1 = (y-1)*iisize + x
               do z=1,nz_fft
                  position = pos1
                  do iy = y,y2
                     do ix=x,x2
                        buf1(position) = buf3(z,ix,iy)
                        position = position+1
                     enddo
                     position = position + n
                  enddo 
                  pos1 = pos1 + iisize*jjsize
               enddo
            enddo
         enddo

         tc = tc + MPI_Wtime()
         
         t = t - MPI_Wtime() 
         call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, &
              buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)

         t = t + MPI_Wtime() 

! Unpack receive buffers into dest

         position=1
         do i=0,jproc-1
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
            enddo
         enddo
         
      endif
#endif

      
      return
      end subroutine
