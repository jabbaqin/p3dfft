! This routine is called only when jproc=1, and only when stride1 is used
! transform backward in Z and transpose array in memory 

      subroutine reorder_trans_b(A,B,C)

      complex(mytype) B(iisize,ny_fft,nz_fft)
      complex(mytype) A(nz_fft,iisize,ny_fft)
      complex(mytype) C(nz_fft,iisize,NB*iisize*nz_fft)
      integer x,y,z,iy,iz,y2,z2

      if(OW) then
         do y=1,ny_fft,nb
            y2 = min(y+nb-1,ny_fft)
            call exec_b_c2(A(1,1,y),1,nz_fft,A(1,1,y),1,nz_fft,nz_fft,nb*iisize)
            do z=1,nz_fft,nb
               z2 = min(z+nb-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = A(iz,x,iy)
                     enddo
                  enddo
               enddo
               
            enddo                  
         enddo

      else
         do y=1,ny_fft,nb
            y2 = min(y+nb-1,ny_fft)
            call exec_b_c2(A(1,1,y),1,nz_fft,C(1,1,y),1,nz_fft,nz_fft,nb*iisize)
            do z=1,nz_fft,nb
               z2 = min(z+nb-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = C(iz,x,iy)
                     enddo
                  enddo
               enddo
               
            enddo                  
         enddo
         
      endif

      return
      end subroutine


! This routine is called only when jproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Z

      subroutine reorder_trans_f(A,B)

      complex(mytype) A(iisize,ny_fft,nz_fft)
      complex(mytype) B(nz_fft,iisize,ny_fft)
!      complex(mytype) C(nz_fft,iisize,ny_fft)
      integer x,y,z,iy,iz,y2,z2

      do y=1,ny_fft,nb
         y2 = min(y+nb-1,ny_fft)
         do z=1,nz_fft,nb
            z2 = min(z+nb-1,nz_fft)
            do iz=z,z2
               do iy=y,y2
                  do x=1,iisize
                     B(iz,x,iy) = A(x,iy,iz)
                  enddo
               enddo
            enddo
         enddo                  
         call exec_f_c2(B(1,1,y),1,nz_fft,B(1,1,y),1,nz_fft,nz_fft,nb*iisize)
      enddo

      return
      end subroutine

