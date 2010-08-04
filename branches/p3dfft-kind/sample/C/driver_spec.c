/*
! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with a 3D sine wave, then 
! performs 3D FFT forward transform, and computes power spectrum.
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
*/

#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double static FORTNAME(t1),FORTNAME(t2),FORTNAME(t3),FORTNAME(t4);
/* double t1,t2,t3,t4,tp1; */
#ifndef SINGLE_PREC
    void compute_spectrum(double *B,int *st,int *sz,int *en,int *ng,double *E,int kmax,int root);
#else
    void compute_spectrum(float *B,int *st,int *sz,int *en,int *ng,float *E,int kmax,int root);
#endif

int main(int argc,char **argv)
{
#ifndef SINGLE_PREC
   double *A,*B, *E,*p1,*p;
#else
   float *A,*B, *E,*p1,*p;
#endif
   int i,j,k,x,y,z,nx,ny,nz,proc_id,nproc,dims[2],ndim,nu;
   int istart[3],isize[3],iend[3];
   int fstart[3],fsize[3],fend[3];
   int iproc,jproc,ng[3],kmax,iex,conf,m,n;
   long int Ntot;

#ifndef SINGLE_PREC
   double pi,twopi,sinyz,cdiff,ccdiff,ans;
   double *sinx,*siny,*sinz,factor,r;
#else
   float pi,twopi,sinyz,cdiff,ccdiff,ans;
   float *sinx,*siny,*sinz,factor,r;
#endif

   double rtime1,rtime2,gt1,gt2,gt3,gt4,gtcomm,tcomm;
   FILE *fp;

#ifndef SINGLE_PREC
   void print_all(double *,long int,int,long int),mult_array(double *,long int,double);
#else
   void print_all(float *,long int,int,long int),mult_array(float *,long int,double);
#endif

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

   pi = atan(1.0)*4.0;
   twopi = 2.0*pi;

   gt1=gt2=gt3=gt4=0.0;

 
   if(proc_id == 0) {
     if((fp=fopen("stdin", "r"))==NULL){
        printf("Cannot open file. Setting to default nx=ny=nz=128, ndim=2, n=1.\n");
        nx=ny=nz=128; n=1;
     } else {
        fscanf(fp,"%d %d %d %d %d\n",&nx,&ny,&nz,&ndim,&n);
        fclose(fp);
     }
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n %d proc. dimensions\n",nx,ny,nz,ndim);
#else
     printf("Single precision\n (%d %d %d) grid\n %d proc. dimensions\n",nx,ny,nz,ndim);
#endif
   }
   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   /*    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);*/
   MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);
   
   if(ndim == 1) {
     dims[0] = 1; dims[1] = nproc;
   }
   else if(ndim == 2) {
     fp = fopen("dims","r");
     if(fp != NULL) {
       if(proc_id == 0)
         printf("Reading proc. grid from file dims\n");
       fscanf(fp,"%d %d\n",dims,dims+1);
       fclose(fp);
       if(dims[0]*dims[1] != nproc) 
          dims[1] = nproc / dims[0];
     }
     else {
       if(proc_id == 0) 
          printf("Creating proc. grid with mpi_dims_create\n");
       dims[0]=dims[1]=0;
       MPI_Dims_create(nproc,2,dims);
       if(dims[0] > dims[1]) {
          dims[0] = dims[1];
          dims[1] = nproc/dims[0];
       }
     }
   }

   if(proc_id == 0) 
      printf("Using processor grid %d x %d\n",dims[0],dims[1]);

   /* Initialize P3DFFT */
   p3dfft_setup(dims,nx,ny,nz,1);
   /* Get dimensions for input and output arrays */
   conf = 1;
   p3dfft_get_dims(istart,iend,isize,conf);
   conf = 2;
   p3dfft_get_dims(fstart,fend,fsize,conf);

#ifndef SINGLE_PREC
   sinx = malloc(sizeof(double)*nx);
   siny = malloc(sizeof(double)*ny);
   sinz = malloc(sizeof(double)*nz);
#else
   sinx = malloc(sizeof(float)*nx);
   siny = malloc(sizeof(float)*ny);
   sinz = malloc(sizeof(float)*nz);
#endif

   for(z=0;z < isize[2];z++)
     sinz[z] = sin((z+istart[2]-1)*twopi/nz);
   for(y=0;y < isize[1];y++)
     siny[y] = sin((y+istart[1]-1)*twopi/ny);
   for(x=0;x < isize[0];x++)
     sinx[x] = sin((x+istart[0]-1)*twopi/nx);

   /* Allocate arrays */

#ifndef SINGLE_PREC
   A = (double *) malloc(sizeof(double) * isize[0]*isize[1]*isize[2]);
   B = (double *) malloc(sizeof(double) * fsize[0]*fsize[1]*fsize[2]*2);
#else
   A = (float *) malloc(sizeof(float) * isize[0]*isize[1]*isize[2]);
   B = (float *) malloc(sizeof(float) * fsize[0]*fsize[1]*fsize[2]*2);
#endif

   /* Initialize array A */
   p1 = A;
   for(z=0;z < isize[2];z++)
     for(y=0;y < isize[1];y++) {
       sinyz = siny[y]*sinz[z];
       for(x=0;x < isize[0];x++)
          *p1++ = sinx[x]*sinyz;
     }

   /*
   for(z=0;z < isize[2];z++)
     for(y=0;y < isize[1];y++) 
       for(x=0;x < isize[0];x++) {
          r = rand()/RAND_MAX;
	  *p1++ =  r;
       }
   */

   Ntot = fsize[0]*fsize[1]*fsize[2]*2;
   factor = 1.0 / (nx*ny) / nz;

   MPI_Barrier(MPI_COMM_WORLD);
   /* Compute forward Fourier transform on A, store results in B */
   rtime1 = - MPI_Wtime();
   p3dfft_ftran_r2c(A,B);
   rtime1 = rtime1 + MPI_Wtime();

   /* normalize */
   mult_array(B,Ntot,factor);

#ifndef SINGLE_PREC
   kmax = sqrt((double) (nx*nx + ny*ny + nz * nz)) *0.5 +0.5;
   E = malloc(sizeof(double) * (kmax+1));
#else
   kmax = sqrt((float) (nx*nx + ny*ny + nz * nz)) *0.5 +0.5;
   E = malloc(sizeof(float) * (kmax+1));
#endif

   ng[0] = nx;	  
   ng[1] = ny;
   ng[2] = nz;

   /* Compute power spectrum, store results in E */
   compute_spectrum(B,fstart,fsize,fend,ng,E,kmax,0);

   if(proc_id == 0) {
     printf("Power spectrum:\n");
     for(i=0; i <= kmax;i++)
       printf("(%d.0 %g)",i,E[i]);
     printf("\n");
   }

   p3dfft_clean();
   free(A);
   free(B);
   free(E);
  
   /* Gather timing statistics */
  MPI_Reduce(&rtime1,&rtime2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  /*
  MPI_Reduce(&FORTNAME(t1),&gt1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t2),&gt2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t3),&gt3,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&FORTNAME(t4),&gt4,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  tcomm = FORTNAME(t1)+FORTNAME(t2)+FORTNAME(t3)+FORTNAME(t4);
  MPI_Reduce(&tcomm,&gtcomm,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  gt1 = gt1 / ((double) n);
  gt2 = gt2 / ((double) n);
  gt3 = gt3 / ((double) n);
  gt4 = gt4 / ((double) n);
  gtcomm = gtcomm / ((double) n);
  */

  if(proc_id == 0) {
     printf("Time per loop=%lg\n",rtime2);
     /*
     printf("Total comm: %g",gtcomm);
     printf("t1=%lg, t2=%lg, t3=%lg, t4=%lg\n",gt1,gt2,gt3,gt4);
     */
  }


  MPI_Finalize();

}

#ifndef SINGLE_PREC
void mult_array(double *A,long int nar,double f)
#else
void mult_array(float *A,long int nar,double f)
#endif
{
  long int i;

  for(i=0;i < nar;i++)
    A[i] *= f;
}

#ifndef SINGLE_PREC
void compute_spectrum(double *B,int *st,int *sz,int *en,int *ng,double *E,int kmax,int root)
#else
void compute_spectrum(float *B,int *st,int *sz,int *en,int *ng,float *E,int kmax,int root)
#endif
  {
    int x,y,z,ik,k2,kx,ky,kz;
#ifndef SINGLE_PREC
    double *el,*p;
#else
    float *el,*p;
#endif

#ifndef SINGLE_PREC
    el = malloc(sizeof(double)*(kmax+1));
#else
    el = malloc(sizeof(float)*(kmax+1));
#endif
    p =B;

    for(ik=0;ik <= kmax;ik++)
      el[ik]=E[ik]=0.0;

    /* Compute local spectrum */
    for(z=0;z < sz[2];z++) {
      kz = z +st[2]-1;
      if(kz > ng[2]/2)
	kz = ng[2] - kz;
      
      for(y=0;y < sz[1];y++) {
	ky = y +st[1]-1;
	if(ky > ng[1]/2)
	  ky = ng[1] - ky;

	for(x=0;x < sz[0]*2;x++) {

	  /* Array B is floats storing complex numbers, so two elements 
	     share one kx */
	  kx = x/2 + (st[0]-1);
    
	  k2 = kx *kx +ky *ky + kz*kz;
#ifndef SINGLE_PREC
	  ik = sqrt((double) k2)+0.5;
#else
          ik = sqrt((float) k2)+0.5;
#endif
	  el[ik] += k2 * (*p) * (*p);
	  p++;
	}
      }
    }

    /* Add up results among processors */
    MPI_Reduce(el,E,kmax+1,MPI_FLOAT,MPI_SUM,root,MPI_COMM_WORLD);

    free(el);

  }

  