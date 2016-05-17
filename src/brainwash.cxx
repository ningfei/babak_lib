///////////////////////////////////////////////////////////////// 
// brainwash.cxx 
// Copyright 2016 by Babak A. Ardekani 
// ALL RIGHTS RESERVED. No part of this program may be 
// used, transferred, or modified by any means without 
// prior written permission. Making copies of any part
// of this program for any purpose is a violation of
// copyright laws.
///////////////////////////////////////////////////////////////// 

#include <stdio.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <dirent.h>
#include "../include/volume.h"
#include "../include/spm_analyze.h"
#include "../include/babak_lib.h"
#include "../include/sph.h"
#include "../include/smooth.h"
#include "../include/minmax.h"
#include "../include/interpolator.h"
#include "../include/stats.h"

#define YES 1
#define NO 0
#define MAXITER 10
#define MAXR 15
#define MAXNATLAS 150

#define MATRIXSIZE 256
#define VOXELSIZE 1.0

//////////////////////////////////////////////////////////////////////////////////////////////////

float FWHM=1.0;
char atlasfilename[MAXNATLAS][1024];

// multi-resolution image dimensions
DIM dim1, dim2, dim4, dim8;

int opt;

static struct option options[] =
{
   {"-version", 0, 'V'},
   {"-Version", 0, 'V'},

   {"-iter8", 1, '8'},
   {"-iter4", 1, '4'},
   {"-iter2", 1, '2'},
   {"-iter1", 1, '1'},

   {"-r", 1, 'r'},
   {"-R", 1, 'R'},

   {"-thresh", 1, 't'},
   {"-threshold", 1, 't'},
   {"-t", 1, 't'},

   {"-sub", 1, 's'},
   {"-sublm", 1, 'M'},

   {"-atlaslist", 1, 'l'},
   {"-list", 1, 'l'},
   {"-l", 1, 'l'},

   {"-atlasdir", 1, 'd'},
   {"-dir", 1, 'd'},
   {"-d", 1, 'd'},

   {"-ppm", 0, 'p'},
   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h",0,'h'},
   {"-help",0,'h'},
   {0, 0, 0}
};

int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

extern float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[]);

void print_help_and_exit()
{
   printf("\n\nUsage:\n"
   "\tbrainwash [-v or -verbose] [-h or -help] [-version]\n"
   "\t[-r <patch radius>] [-R <search radius>] [-t <percent threshold>]\n"
   "\t-sub <subject image>\n\n" 

   "Required arguments:\n\n"
   "\t-sub <subject image>\n"
   "\t\tSpecifies the image to be registred to the <target image>. <subject image> is expected\n"
   "\t\tto be of type short in NIFTI format.\n\n"

   "Optional arguments:\n\n"
   "\t-v or -verbose\n"
   "\t\tEnables verbose mode.\n\n"

   "\t-r <patch radius>\n"
   "\t\tSpecifies the patch radius (default=3 pixels).\n\n"

   "\t-R <search radius>\n"
   "\t\tSpecifies the search radius (default=3 pixels).\n\n"

   "\t-h or -help\n"
   "\t\tPrints help message.\n\n"

   "\t-version\n"
   "\t\tPrints program vesion.\n\n"
   );
   exit(0);
}

short *computeReslicedImage2(short *im1, DIM dim1, DIM dim2, float *Xwarp, float *Ywarp, float *Zwarp)
{
 	float  x,y,z;   
	int q;
	int np1;
	short *im2;
	float xc1, yc1, zc1;
	float xc2, yc2, zc2;

	float *beta, del;
	float *c;

	if(opt_cubicspline)
	{
		beta=computeBeta(&del);
		c = (float *)calloc(dim1.nv, sizeof(float));
		cubicSplineAnalysis(im1, c, dim1.nx, dim1.ny, dim1.nz);
	}

	im2=(short *)calloc(dim2.nv,sizeof(short));

	xc1=dim1.dx*(dim1.nx-1)/2.0;     /* +---+---+ */
	yc1=dim1.dy*(dim1.ny-1)/2.0;
	zc1=dim1.dz*(dim1.nz-1)/2.0;

   	xc2=dim2.dx*(dim2.nx-1)/2.0;     /* +---+---+ */
	yc2=dim2.dy*(dim2.ny-1)/2.0;
	zc2=dim2.dz*(dim2.nz-1)/2.0;

	q=0;
	for(int k=0;k<dim2.nz;k++) 
	for(int j=0;j<dim2.ny;j++) 
  	for(int i=0;i<dim2.nx;i++) 
	{
		z = (k*dim2.dz - zc2 + Zwarp[q] + zc1) /dim1.dz;
		y = (j*dim2.dy - yc2 + Ywarp[q] + yc1) /dim1.dy;
		x = (i*dim2.dx - xc2 + Xwarp[q] + xc1) /dim1.dx;

		if(opt_cubicspline)
			im2[q++] = (short)(cubicSplineSynthesis(c, dim1.nx, dim1.ny, dim1.nz, x, y, z, beta, del)+0.5);
		else
			im2[q++]=(short)(linearInterpolator(x,y,z,im1,dim1.nx,dim1.ny,dim1.nz,dim1.np)+0.5);
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

void brainwashnlReg(short *sub, short *trg, short *msk, DIM dim, int r, int R, float *Xwarp, float *Ywarp, float *Zwarp, int S)
{
   float I[16];
   float *Xw, *Yw, *Zw;

   SPH subsph(r);
   SPH trgsph(r);
   SPH searchsph(R);

   SHORTIM subim;
   set_dim(subim,dim);
   subim.v=sub;

   SHORTIM trgim;
   set_dim(trgim,dim);
   trgim.v=trg;

   int v; // voxel index going from 0 to dim.nv-1
   int P[3], Q[3];

   Xw = (float *)calloc(dim.nv, sizeof(float));
   Yw = (float *)calloc(dim.nv, sizeof(float));
   Zw = (float *)calloc(dim.nv, sizeof(float));

   for(int i=0; i<dim.nv; i++) if(msk[i]==0) sub[i]=trg[i]=0;
    
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np + j*dim.nx + i;

      if(msk[v]==100 || msk[v]<=0 || trg[v]<=0)
      {
         continue;
      }

      trgsph.set(trgim,i,j,k);
      standardize(trgsph.v,trgsph.v,trgsph.n);

      P[0]=i; P[1]=j; P[2]=k;
      detect_lm(searchsph, subsph, subim, P, trgsph, Q);

      // make sure Q is within the mask
      if(msk[ Q[2]*dim.np + Q[1]*dim.nx + Q[0] ] != 0)
      {
         Xw[v] = (Q[0]-P[0])*dim.dx;
         Yw[v] = (Q[1]-P[1])*dim.dy;
         Zw[v] = (Q[2]-P[2])*dim.dz;
      }
   }

   float *Xww, *Yww, *Zww;
   Xww = (float *)calloc(dim.nv, sizeof(float));
   Yww = (float *)calloc(dim.nv, sizeof(float));
   Zww = (float *)calloc(dim.nv, sizeof(float));

   {
      SPH xsph(S);
      SPH ysph(S);
      SPH zsph(S);

      for(int k=0; k<dim.nz; k++)
      for(int j=0; j<dim.ny; j++)
      for(int i=0; i<dim.nx; i++)
      {
         v = k*dim.np + j*dim.nx + i;

         if(trg[v]<=0)
         {
            continue;
         }

         xsph.set(Xw,dim.nx,dim.ny,dim.nz,i,j,k);
         ysph.set(Yw,dim.nx,dim.ny,dim.nz,i,j,k);
         zsph.set(Zw,dim.nx,dim.ny,dim.nz,i,j,k);

         Xww[v] = xsph.mean();
         Yww[v] = ysph.mean();
         Zww[v] = zsph.mean();
      }

      // replace Xw with Xww etc.
      free(Xw); free(Yw); free(Zw);
      Xw=Xww; Yw=Yww; Zw=Zww;
   }

   float *tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Xw, dim, dim1, I, LIN); 
   free(Xw); Xw=tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Yw, dim, dim1, I, LIN); 
   free(Yw); Yw=tmp;

   set_to_I(I,4); 
   tmp = resliceImage(Zw, dim, dim1, I, LIN); 
   free(Zw); Zw=tmp;

   for(int i=0; i<dim1.nv; i++) {
      Xwarp[i] += Xw[i];
      Ywarp[i] += Yw[i];
      Zwarp[i] += Zw[i];
   }
   free(Xw); free(Yw); free(Zw);
}

// finds an affine registration in "A" to takes points from sub to trg space
void affineReg(short *sub, short *trg, short *msk, DIM dim, int r, int R, float *A)
{
   SPH subsph(r);
   SPH trgsph(r);
   SPH searchsph(R);

   SHORTIM subim;
   set_dim(subim,dim);
   subim.v=sub;

   SHORTIM trgim;
   set_dim(trgim,dim);
   trgim.v=trg;

   int *PT, *QT;
   char *flg;

   float *P, *Q;

   flg = (char *)calloc(dim.nv, sizeof(char));
   PT = (int *)calloc(dim.nv*3, sizeof(int));
   QT = (int *)calloc(dim.nv*3, sizeof(int));
   
   int c=0;
   int v;
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np+ j*dim.nx+ i;

      if(msk[v]<=0 || sub[v]<=0 || trg[v]<=0)
      {
         flg[v]=0;
         continue;
      }

      subsph.set(subim,i,j,k);
      standardize(subsph.v,subsph.n);

      PT[3*c + 0] = QT[3*c + 0] = i;
      PT[3*c + 1] = QT[3*c + 1] = j;
      PT[3*c + 2] = QT[3*c + 2] = k;
      detect_lm(searchsph, trgsph, trgim, PT+3*c, subsph, QT+3*c);

      flg[v]=1;
      c++;
   }

   int n=c;

   P = (float *)calloc(3*n, sizeof(float));
   Q = (float *)calloc(3*n, sizeof(float));

   c=0;
   for(int i=0; i<dim.nv; i++)
   if(flg[i])
   {
      P[0*n + c]=(PT[3*c + 0] - (dim.nx-1)/2.0)*dim.dx;
      P[1*n + c]=(PT[3*c + 1] - (dim.ny-1)/2.0)*dim.dy;
      P[2*n + c]=(PT[3*c + 2] - (dim.nz-1)/2.0)*dim.dz;

      Q[0*n + c]=(QT[3*c + 0] - (dim.nx-1)/2.0)*dim.dx;
      Q[1*n + c]=(QT[3*c + 1] - (dim.ny-1)/2.0)*dim.dy;
      Q[2*n + c]=(QT[3*c + 2] - (dim.nz-1)/2.0)*dim.dz;

      c++;
   }

   leastSquaresAffineTrans(P, Q, n, A);

   free(PT); 
   free(QT);
   free(flg);
   free(P); 
   free(Q);
}

void generateMultiResolution(short *sub, DIM subdim, float *T, float *Xwarp, float *Ywarp, float *Zwarp, short *im1, short *im2, short *im4, short *im8)
{
   short *tmp;

   float I[16];
   float *Xw,*Yw,*Zw;

   Xw = (float *)calloc(dim1.nv, sizeof(float));
   Yw = (float *)calloc(dim1.nv, sizeof(float));
   Zw = (float *)calloc(dim1.nv, sizeof(float));

   for(int i=0; i<dim1.nv; i++) {
      Xw[i] = Xwarp[i];
      Yw[i] = Ywarp[i];
      Zw[i] = Zwarp[i];
   }

   combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xw, Yw, Zw, T);

   tmp=computeReslicedImage2(sub, subdim, dim1, Xw, Yw, Zw);
   for(int i=0; i<dim1.nv; i++) im1[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im1, dim1, dim2, I, LIN); 
   for(int i=0; i<dim2.nv; i++) im2[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im2, dim2, dim4, I, LIN); 
   for(int i=0; i<dim4.nv; i++) im4[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im4, dim4, dim8, I, LIN); 
   for(int i=0; i<dim8.nv; i++) im8[i]=tmp[i];
   free(tmp);

   free(Xw); free(Yw); free(Zw);
}

void generateMultiResolution(short *im, DIM im_dim, float *T, short *im1, short *im2, short *im4, short *im8)
{
   short *tmp;
   float I[16];
   float *invT;		

   invT=inv4(T); 
   tmp = resliceImage(im, im_dim, dim1, invT, LIN); 
   free(invT);
   for(int i=0; i<dim1.nv; i++) im1[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im1, dim1, dim2, I, LIN); 
   for(int i=0; i<dim2.nv; i++) im2[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im2, dim2, dim4, I, LIN); 
   for(int i=0; i<dim4.nv; i++) im4[i]=tmp[i];
   free(tmp);

   set_to_I(I,4); 
   tmp = resliceImage(im4, dim4, dim8, I, LIN); 
   for(int i=0; i<dim8.nv; i++) im8[i]=tmp[i];
   free(tmp);

   return;
}

void combine_warps_and_trans(float *trgTPIL, float *Xwarp, float *Ywarp, float *Zwarp, float *Xout, float *Yout, float *Zout, DIM trgdim)
{
   int v;
   float xc,yc,zc;
   float xc1,yc1,zc1;
   float x,y,z;   
   float x1,y1,z1;   
   float i1,j1,k1;   

   xc=trgdim.dx*(trgdim.nx-1.0)/2.0;     /* +---+---+ */
   yc=trgdim.dy*(trgdim.ny-1.0)/2.0;
   zc=trgdim.dz*(trgdim.nz-1.0)/2.0;

   xc1=dim1.dx*(dim1.nx-1.0)/2.0;     /* +---+---+ */
   yc1=dim1.dy*(dim1.ny-1.0)/2.0;
   zc1=dim1.dz*(dim1.nz-1.0)/2.0;
 
   for(int k=0;k<trgdim.nz;k++) 
   for(int j=0;j<trgdim.ny;j++) 
   for(int i=0;i<trgdim.nx;i++) 
   {
      v = k*trgdim.np + j*trgdim.nx + i;

      // (i*dx-xc) converts from image coordinates (i,j,z) to (x,y,z) coordinates
      x = (i*trgdim.dx - xc);
      y = (j*trgdim.dy - yc);
      z = (k*trgdim.dz - zc);

      x1 =  trgTPIL[0]*x +trgTPIL[1]*y +trgTPIL[2]*z  +trgTPIL[3];
      y1 =  trgTPIL[4]*x +trgTPIL[5]*y +trgTPIL[6]*z  +trgTPIL[7];
      z1 =  trgTPIL[8]*x +trgTPIL[9]*y +trgTPIL[10]*z +trgTPIL[11];

      i1 = (x1 + xc1) / dim1.dx;
      j1 = (y1 + yc1) / dim1.dy;
      k1 = (z1 + zc1) / dim1.dz;

      x1 += linearInterpolator(i1,j1,k1,Xwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      y1 += linearInterpolator(i1,j1,k1,Ywarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      z1 += linearInterpolator(i1,j1,k1,Zwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      
      Xout[v] = x1 - i*trgdim.dx + xc;
      Yout[v] = y1 - j*trgdim.dy + yc;
      Zout[v] = z1 - k*trgdim.dz + zc;
   }
}

void approximate_affine(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
	int np;
	int q,N;
	float xc, yc, zc;
	float rx,ry,rz;
	float sx,sy,sz;
  	float  x,y,z;   
	float *AFF, *invAFF;	// affine transform

	double Mrx, Mry, Mrz;
	double Msx, Msy, Msz;
	double SR[9], RR[9];
	double *invRR;
	double A[9],B[3];

   AFF = (float *)calloc(16,sizeof(float));
   if(AFF==NULL)
   {
      printf("\nMemory allocation error (AFF) ...\n");
      exit(1);
   }

	np = nx*ny;

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	Mrx=Mry=Mrz=0.0;
	Msx=Msy=Msz=0.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	N = 0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			Mrx += rx; Mry += ry; Mrz += rz;
			Msx += sx; Msy += sy; Msz += sz;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N; Mrz /= N;
		Msx /= N; Msy /= N; Msz /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 3x3 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<9; i++) SR[i]=RR[i]=0.0;

	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			rx -= Mrx; ry -= Mry; rz -= Mrz;
			sx -= Msx; sy -= Msy; sz -= Msz;

			SR[0]+=sx*rx; SR[1]+=sx*ry; SR[2]+=sx*rz;
			SR[3]+=sy*rx; SR[4]+=sy*ry; SR[5]+=sy*rz;
			SR[6]+=sz*rx; SR[7]+=sz*ry; SR[8]+=sz*rz;

			RR[0]+=rx*rx; RR[1]+=rx*ry; RR[2]+=rx*rz;
			RR[3]+=ry*rx; RR[4]+=ry*ry; RR[5]+=ry*rz;
			RR[6]+=rz*rx; RR[7]+=rz*ry; RR[8]+=rz*rz;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv3(RR);
	multi(SR,3,3,invRR,3,3,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry - A[2]*Mrz;
	B[1] = Msy - A[3]*Mrx - A[4]*Mry - A[5]*Mrz;
	B[2] = Msz - A[6]*Mrx - A[7]*Mry - A[8]*Mrz;

	// Eq. (3) of tech. notes
	AFF[0]=(float)A[0]; AFF[1]=(float)A[1];  AFF[2]=(float)A[2];  AFF[3]=(float)B[0];
	AFF[4]=(float)A[3]; AFF[5]=(float)A[4];  AFF[6]=(float)A[5];  AFF[7]=(float)B[1];
	AFF[8]=(float)A[6]; AFF[9]=(float)A[7]; AFF[10]=(float)A[8]; AFF[11]=(float)B[2];
	AFF[12]=0.0; AFF[13]=0.0; AFF[14]=0.0; AFF[15]=1.0;

	// Eq. (3.5) of tech. notes
	invAFF = inv4(AFF);
	delete AFF;

	///////////////////////////////////////////////////////////////////////////////

	// replace T with invAFF
	for(int i=0; i<16; i++) T[i] = invAFF[i];

	delete invAFF;
}

void fillzerovectors(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
  	float  x,y,z;   
  	float  xx,yy,zz;   
	float xc,yc,zc;
	float *invT;		
	int q, np;

    np = nx*ny;

	invT=inv4(T);

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	q=0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]==0.0 && Ywarp[q]==0.0 && Zwarp[q]==0.0)
		{
		   xx = (i*dx - xc);
		   yy = (j*dy - yc);
		   zz = (k*dz - zc);

		   x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]  );
		   y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]  );
		   z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11] );

		   Xwarp[q] = x - xx;
		   Ywarp[q] = y - yy;
		   Zwarp[q] = z - zz;
	   }
	}

	free(invT);
}

void ivf(float *Xwarp1, float *Ywarp1, float *Zwarp1, DIM dim1, float *Xwarp2, float *Ywarp2, float *Zwarp2, DIM dim2)
{
   int i2, j2, k2;
   int ii, jj, kk;
   int np1, nv1;
   int np2, nv2;
   int wx, wy, wz;
   int v1, v1_part1, v1_part2;
   int v2, v2_part1, v2_part2;

   float sigma, K;
   float *W, w;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float x,y,z;
   float xx, yy, zz;
   float x2, y2, z2;
   float x1, y1, z1;

   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   ////////////////////////////////////////////////////////////////////////////////

   nx1 = dim1.nx; ny1 = dim1.ny; nz1 = dim1.nz;
   nx2 = dim2.nx; ny2 = dim2.ny; nz2 = dim2.nz;

   dx1 = dim1.dx; dy1 = dim1.dy; dz1 = dim1.dz;
   dx2 = dim2.dx; dy2 = dim2.dy; dz2 = dim2.dz;

   ////////////////////////////////////////////////////////////////////////////////

   nv1 = nx1 * ny1 * nz1;
   np1 = nx1 * ny1;

   nv2 = nx2 * ny2 * nz2;
   np2 = nx2 * ny2;

   ////////////////////////////////////////////////////////////////////////////////

   sigma = FWHM/2.35482; // converts FWHM to standard deviation
   K = 2.0*sigma*sigma;

   // Twice the standard deviation is greater than 95% of the area
   wx = (int)ceilf(2.0*sigma/dx2);
   wy = (int)ceilf(2.0*sigma/dy2);
   wz = (int)ceilf(2.0*sigma/dz2);

   ////////////////////////////////////////////////////////////////////////////////

//   if(opt_v)
//   {
//      printf("Standard deviation of the Gaussian kernel = %f\n\n",sigma);
//      printf("Support region size: %d x %d x %d voxels\n\n", 2*wx+1, 2*wy+1, 2*wz+1);

//      printf("Input matrix size = %d x %d x %d\n",nx1, ny1, nz1);
//      printf("Input voxel size = %7.5f x %7.5f x %7.5f mm^3\n\n",dx1, dy1, dz1);

//      printf("Output matrix size = %d x %d x %d\n",nx2, ny2, nz2);
//      printf("Output voxel size = %7.5f x %7.5f x %7.5f mm^3\n",dx2, dy2, dz2);
//   }

   ////////////////////////////////////////////////////////////////////////////////

   W = (float *)calloc(nv2,sizeof(float));
   if(W==NULL)
   {
      printf("\nMemory allocation error (W) ...\n");
      exit(1);
   }

   xc1 = dx1 * (nx1-1.0)/2.0;
   yc1 = dy1 * (ny1-1.0)/2.0;
   zc1 = dz1 * (nz1-1.0)/2.0;

   xc2 = dx2 * (nx2-1.0)/2.0;
   yc2 = dy2 * (ny2-1.0)/2.0;
   zc2 = dz2 * (nz2-1.0)/2.0;

   for(int k1=0; k1<nz1; k1++)
   {
      v1_part1 = k1*np1;
      z1 = dz1 * k1 - zc1;
      
      for(int j1=0; j1<ny1; j1++)
      {
         v1_part2 = j1*nx1;
         y1 = dy1 * j1 - yc1;

         for(int i1=0; i1<nx1; i1++)
         {
            v1 = v1_part1 + v1_part2 + i1;
            x1 = dx1 * i1 - xc1;

            if( Xwarp1[v1]!=0.0 && Ywarp1[v1]!=0.0 && Zwarp1[v1]!=0.0 )
            {
               x2 = x1 + Xwarp1[v1];
               y2 = y1 + Ywarp1[v1];
               z2 = z1 + Zwarp1[v1];

               i2 = (int)nearbyintf( (x2 + xc2)/dx2 );
               j2 = (int)nearbyintf( (y2 + yc2)/dy2 );
               k2 = (int)nearbyintf( (z2 + zc2)/dz2 );

               for(int k=-wz; k<=wz; k++)
               {
                  kk = k+k2;
                  v2_part1 = kk*np2;

                  for(int j=-wy; j<=wy; j++)
                  {
                     jj = j+j2;
                     v2_part2 = jj*nx2;

                     for(int i=-wx; i<=wx; i++)
                     {
                        ii = i+i2;

                        if( ii>0 && jj>0 && kk>0 && ii<nx2 && jj<ny2 && kk<nz2 )
                        {
                           x = dx2 * ii - xc2;
                           y = dy2 * jj - yc2;
                           z = dz2 * kk - zc2;

                           v2 = v2_part1 + v2_part2 + ii;

                           xx = x-x2;
                           yy = y-y2;
                           zz = z-z2;

                           w = expf( -(xx*xx + yy*yy + zz*zz)/K);

                           W[v2] += w;

                           Xwarp2[v2] += w*( x1 - x);
                           Ywarp2[v2] += w*( y1 - y);
                           Zwarp2[v2] += w*( z1 - z);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   for(int i=0; i<nv2; i++)
   if(W[i]>0.0)
   {
      Xwarp2[i] /= W[i];
      Ywarp2[i] /= W[i];
      Zwarp2[i] /= W[i];
   }

   {
      float T[16];
      approximate_affine(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
      fillzerovectors(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
   }

   free(W);
}

void read_default_atlas_names(const char *brainwashatlasdir, int &natlas)
{
   char fileprefix[1024];

   DIR *dp;
   // struct direct *dir; // this one didn't work on SUN 
   struct dirent *dir; 

   natlas=0;

   // open the directory and associate the directory stream dp with it 
   dp=opendir(brainwashatlasdir);

   // If opendir() returns a NULL pointer, then brainwashatlasdir cannot be accessed, 
   // or if it cannot malloc() enough memory to hold the whole thing. 
   if(dp==NULL)
   {
      printf("Warning: Cannot access %s, returning ...\n",brainwashatlasdir);
      return;
   }

   // Get a pointer to the next directory entry using readdir().  If readdir()
   // returns NULL, the end of the directory has been reached.  
   while( (dir=readdir(dp)) != NULL && natlas<MAXNATLAS)
   {
      // If the string length is less than 5, it cannot be a NIFTI file since
      // NIFTI files have names like *.nii which are at least 5 characters
      if(strlen(dir->d_name) < 5 )
         continue;

      // then looks like a real NIFTI image since niftiFilename does some checks
      if( niftiFilename(fileprefix, dir->d_name)==1 ) 
      {
         sprintf(atlasfilename[natlas],"%s/%s",brainwashatlasdir,dir->d_name);
         natlas++; // increment natlas to indicate that one more atlas has been read 
      }
   };

   // close the directory stream and free the structure associated with dp 
   closedir(dp);
}

void read_atlas_list(const char *atlaslistfile, int &natlas)
{
   FILE *fp;
   char filename[1024];
   char fileprefix[1024];

   natlas=0;

   fp=fopen(atlaslistfile,"r");

   if(fp==NULL)
   {
      printf("Warning: Cannot access %s, returning ...\n",atlaslistfile);
      return;
   }

   while( fscanf(fp,"%s",filename) != EOF && natlas<MAXNATLAS)
   {
      if(strlen(filename) < 5 )
         continue;

      // then looks like a real NIFTI image since niftiFilename does some checks
      if( niftiFilename(fileprefix, filename)==1 ) 
      {
         sprintf(atlasfilename[natlas],"%s",filename);
         natlas++; // increment natlas to indicate that one more atlas has been read 
      }
   };

   fclose(fp);
}

int main(int argc, char **argv)
{
   float threshold=50.0;

   char atlaslistfile[1024]=""; 

   float *label;

   char brainwashatlasdir[1024]=""; 
   int natlas; // number of available atlases
   
   // input "subject image" full path
   char subjectImageFile[1024]=""; 

   // points to the current atlas
   char *targetImageFile;

   // extracted subject image filename only without suffix
   char subprefix[1024]=""; 
   // extracted target image filename only without suffix
   char trgprefix[1024]="";

   // original subject and target images loaded 
   // from subjectImageFile and targetImageFile
   short *sub, *trg;

   // NIFTI headers for input subject and target images
   nifti_1_header sub_hdr, trg_hdr;

   // structures for holding subject and target image dimensions
   DIM subdim, trgdim;

   // a linear transformation that brings the subject/target image from its native space to PIL space
   float subTPIL[16], trgTPIL[16];

   // filenames where AC, PC, and VSPS are manually specific for subject/target image
   char sublmfile[1024]="";
   char trglmfile[1024]=""; // always going to be null

   // number if iterations of NL registration at each resolution level
   int iter8=4;
   int iter4=3;
   int iter2=2;
   int iter1=0;

   // Subject and target PIL versions and brain mask at different resolutons
   short *subPIL1=NULL, *subPIL2=NULL, *subPIL4=NULL, *subPIL8=NULL;
   short *trgPIL1=NULL, *trgPIL2=NULL, *trgPIL4=NULL, *trgPIL8=NULL;
   short *mskPIL1=NULL, *mskPIL2=NULL, *mskPIL4=NULL, *mskPIL8=NULL;

   // intermediate displacement field from trgPIL to subPIL
   float *Xwarp=NULL, *Ywarp=NULL, *Zwarp=NULL;
   float *invXout, *invYout, *invZout;

   // 4x4 identity matrix
   float I[16];

   // space holder for a generic filename
   char filename[1024];

   int patch_r=3; // patch radius
   int search_R=3; // search radius 

   // subject to target image affine transformation
   float sub_to_trg[16];

   opt_ppm = NO;
   opt_txt=NO;

   /////////////////////////////////////////////////////////////////////
   // maximize the stack size
   /////////////////////////////////////////////////////////////////////
   {
      // define a rlimit structure
      struct rlimit rlpty;

      // get the resource limit of type RLIMIT_STACK for stack size
      (void) getrlimit(RLIMIT_STACK,&rlpty);

      // set the stack size to its maximum
      rlpty.rlim_cur=rlpty.rlim_max;
      (void) setrlimit(RLIMIT_STACK,&rlpty);
   }
   /////////////////////////////////////////////////////////////////////

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'V':
            printf("Version 2.0 (last modified: May 17, 2016)\n");
            exit(0);
         case 't':
            threshold = atof(optarg);
            break;
         case 'r':
            patch_r = atoi(optarg);
            break;
         case 'p':
            opt_ppm=YES;
            break;
         case 'R':
            search_R = atoi(optarg);
            break;
         case '8':
            iter8 = atoi(optarg);
            break;
         case '4':
            iter4 = atoi(optarg);
            break;
         case '2':
            iter2 = atoi(optarg);
            break;
         case '1':
            iter1 = atoi(optarg);
            break;
         case 'l':
            sprintf(atlaslistfile,"%s",optarg);
            break;
         case 'd':
            sprintf(brainwashatlasdir,"%s",optarg);
            break;
         case 'M':
            sprintf(sublmfile,"%s",optarg);
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'v':
            opt_v=YES;
            break;
         case 's':
            sprintf(subjectImageFile,"%s",optarg);
            break;
         case '?':
            print_help_and_exit();
		}
	}

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Variable initializations

   getARTHOME();

   // revert to the default value if unreasonable input
   if(threshold<0.0 || threshold>100.0) threshold=50.0;

   // Initialize dim1, dim2, dim4, and dim8 multi-resolution image dimensions
   dim1.nx = dim1.ny = dim1.nz = MATRIXSIZE;
   dim1.np=dim1.nx*dim1.ny; 
   dim1.nv=dim1.np*dim1.nz; 
   dim1.dx = dim1.dy = dim1.dz = VOXELSIZE;

   dim2.nx = dim2.ny = dim2.nz = MATRIXSIZE/2;
   dim2.np=dim2.nx*dim2.ny; 
   dim2.nv=dim2.np*dim2.nz; 
   dim2.dx = dim2.dy = dim2.dz = VOXELSIZE*2.0;

   dim4.nx = dim4.ny = dim4.nz = MATRIXSIZE/4;
   dim4.np=dim4.nx*dim4.ny; 
   dim4.nv=dim4.np*dim4.nz; 
   dim4.dx = dim4.dy = dim4.dz = VOXELSIZE*4.0;

   dim8.nx = dim8.ny = dim8.nz = MATRIXSIZE/8;
   dim8.np=dim8.nx*dim8.nx; 
   dim8.nv=dim8.np*dim8.nz; 
   dim8.dx = dim8.dy = dim8.dz = VOXELSIZE*8.0;

   trgPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   subPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   mskPIL1 = (short *)calloc(dim1.nv,sizeof(short));
   trgPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   subPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   mskPIL2 = (short *)calloc(dim2.nv,sizeof(short));
   trgPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   subPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   mskPIL4 = (short *)calloc(dim4.nv,sizeof(short));
   trgPIL8 = (short *)calloc(dim8.nv,sizeof(short));
   subPIL8 = (short *)calloc(dim8.nv,sizeof(short));
   mskPIL8 = (short *)calloc(dim8.nv,sizeof(short));

   Xwarp = (float *)calloc(dim1.nv, sizeof(float));
   Ywarp = (float *)calloc(dim1.nv, sizeof(float));
   Zwarp = (float *)calloc(dim1.nv, sizeof(float));

   // reset unreasonable values to default
   if(iter8 < 0 || iter8>MAXITER ) iter8=4;
   if(iter4 < 0 || iter4>MAXITER ) iter4=3;
   if(iter2 < 0 || iter2>MAXITER ) iter2=2;
   if(iter1 < 0 || iter1>MAXITER ) iter1=0;

   if(patch_r <= 0 || patch_r > MAXR ) patch_r=3;
   if(search_R <= 0 || search_R > MAXR ) search_R=3;

   if(brainwashatlasdir[0]=='\0')
   {
      sprintf(brainwashatlasdir,"%s/brainwashatlas",ARTHOME);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////

   if(atlaslistfile[0]=='\0')
   {
      if(opt_v) printf("Brainwash atlases are read from %s directory.\n",brainwashatlasdir);
      read_default_atlas_names(brainwashatlasdir, natlas);
   }
   else
   {
      if(opt_v) printf("Brainwash atlases are read from %s list.\n",atlaslistfile);
      read_atlas_list(atlaslistfile, natlas);
   }

   if(opt_v) 
   {
      printf("Number of available atlases = %d \n",natlas);
   }

   if(natlas<=0)
   {
      exit(0);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////

   // compute mskPIL1, mskPIL2, mskPIL4 and mskPIL8
   {
      int2 *PILbraincloud;
      DIM PILbraincloud_dim;
      nifti_1_header PILbraincloud_hdr;

      sprintf(filename,"%s/PILbrain.nii",ARTHOME);

      PILbraincloud = (int2 *)read_nifti_image(filename, &PILbraincloud_hdr);

      if(PILbraincloud==NULL)
      {
            printf("Error reading %s, aborting ...\n", filename);
            exit(1);
      }

      set_dim(PILbraincloud_dim, PILbraincloud_hdr);

      set_to_I(I,4);
      generateMultiResolution(PILbraincloud, PILbraincloud_dim, I, mskPIL1, mskPIL2, mskPIL4, mskPIL8);

      free(PILbraincloud);
   }

   // ensure the subject image is specified at the command line
   if( subjectImageFile[0] == '\0')
   {
      printf("\nPlease specify a \"subject image\" using -sub <filename.nii> ...\n\n");		
      exit(0);
   } 

   // extract the subject filename without path/suffix
   if( niftiFilename(subprefix, subjectImageFile)==0 ) { exit(0); }

   // read the subject image
   sub = (int2 *)read_nifti_image(subjectImageFile, &sub_hdr);
   if(sub==NULL) 
   {
      printf("\nError: Reading subject image %s failed.\n\n",subjectImageFile);
      exit(0);
   }
   set_dim(subdim,sub_hdr); // transfer info from sub_hdr to subdim

   // print some info about the subject image
   if(opt_v)
   {
      printf("Subject image file: %s\n",subjectImageFile);
      printf("\tMatrix size: %d x %d x %d\n",subdim.nx, subdim.ny, subdim.nz);
      printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",subdim.dx, subdim.dy, subdim.dz);
   }

   // find subTPIL using automated MSP, AC/PC and 8 MSP landmarks detection
   if(opt_v) printf("Computing subject image PIL transformation ...\n");
   if(sublmfile[0] != '\0' && opt_v) printf("Subject image landmarks are read from %s\n",sublmfile);
   new_PIL_transform(subjectImageFile, sublmfile, subTPIL);

   if(opt_v) printf("Patch radius = %d mm\n", patch_r);
   if(opt_v) printf("Search radius = %d mm\n", search_R);

   invXout = (float *)calloc(subdim.nv,sizeof(float));
   invYout = (float *)calloc(subdim.nv,sizeof(float));
   invZout = (float *)calloc(subdim.nv,sizeof(float));

   label = (float *)calloc(subdim.nv,sizeof(float));

   for(int a=0; a<natlas; a++)
   {
      float *invT;		
      float A[16];
      short *tmp;

      // final displacement field from trg to sub
      float *Xout=NULL, *Yout=NULL, *Zout=NULL; 

      targetImageFile=atlasfilename[a];

      // extract the target filename without path/suffix
      if( niftiFilename(trgprefix, targetImageFile)==0 ) { exit(0); }

      // read the target image
      trg = (int2 *)read_nifti_image(targetImageFile, &trg_hdr);
      if(trg==NULL) 
      {
         printf("\nError: Reading target image %s failed.\n\n",targetImageFile);
         exit(0);
      }
      set_dim(trgdim,trg_hdr); // transfer info from trg_hdr to trgdim

      // print some info about the target image
      if(opt_v)
      {
         printf("Atlas %d image file: %s\n",a+1,targetImageFile);
         printf("\tMatrix size: %d x %d x %d\n",trgdim.nx, trgdim.ny, trgdim.nz);
         printf("\tVoxel size: %6.4f x %6.4f x %6.4f\n",trgdim.dx, trgdim.dy, trgdim.dz);
      }

      // find trgTPIL using automated MSP, AC/PC and 8 MSP landmarks detection
      if(opt_v) printf("Computing atlas image PIL transformation ...\n");
      new_PIL_transform(targetImageFile, trglmfile, trgTPIL);
   
      // compute trgPIL1, trgPIL2, trgPIL4 and trgPIL8
      generateMultiResolution(trg, trgdim, trgTPIL, trgPIL1, trgPIL2, trgPIL4, trgPIL8);

      ////////////////////////////////////////////////////////////////////////////////////////////
      // Compute the affine transformation  sub_to_trg and update subTPIL
      if(opt_v) printf("Affine registration @ 12.5\% resolution ...\n");
      generateMultiResolution(sub, subdim, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL8, trgPIL8, mskPIL8, dim8, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 25\% resolution ...\n");
      generateMultiResolution(sub, subdim, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL4, trgPIL4, mskPIL4, dim4, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 50\% resolution ...\n");
      generateMultiResolution(sub, subdim, subTPIL, subPIL1, subPIL2, subPIL4, subPIL8);
      affineReg(subPIL2, trgPIL2, mskPIL2, dim2, patch_r, search_R, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      invT = inv4(trgTPIL);
      multi(invT, 4,4, subTPIL, 4, 4, sub_to_trg);
      free(invT);

      if(opt_v) printMatrix(sub_to_trg,4,4,"Atlas image -> target image affine transformation",NULL);
      ////////////////////////////////////////////////////////////////////////////////////////////

      for(int i=0; i<dim1.nv; i++) Xwarp[i]=Ywarp[i]=Zwarp[i]=0.0;

      //////////////////////////////////////////////////////////////////
      // Non-linear registration from trgPIL to subPIL
      
      if(iter8>=1)
      {
         if(opt_v) printf("Non-linear registration @ 12.5\% resolution ...\n");
         for(int i=0; i<iter8; i++)
         {
            generateMultiResolution(sub, subdim, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
            brainwashnlReg(subPIL8, trgPIL8, mskPIL8, dim8, patch_r, search_R/3, Xwarp, Ywarp, Zwarp,5);
         }
      }

      if(iter4>=1)
      {
         if(opt_v) printf("Non-linear registration @ 25\% resolution ...\n");
         for(int i=0; i<iter4; i++)
         {
            generateMultiResolution(sub, subdim, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
            brainwashnlReg(subPIL4, trgPIL4, mskPIL4, dim4, patch_r, 2*search_R/3, Xwarp, Ywarp, Zwarp,5);
         }
      }

      if(iter2>=1)
      {
         if(opt_v) printf("Non-linear registration @ 50\% resolution ...\n");
         for(int i=0; i<iter2; i++)
         {
            generateMultiResolution(sub, subdim, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
            brainwashnlReg(subPIL2, trgPIL2, mskPIL2, dim2, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
         }
      }

      // default value of iter1=0 for brainwash
      if(iter1>=1)
      {
         if(opt_v) printf("Non-linear registration @ 100\% resolution ...\n");
         for(int i=0; i<iter1; i++)
         {
            generateMultiResolution(sub, subdim, subTPIL, Xwarp, Ywarp, Zwarp, subPIL1, subPIL2, subPIL4, subPIL8);
            brainwashnlReg(subPIL1, trgPIL1, mskPIL1, dim1, patch_r, search_R, Xwarp, Ywarp, Zwarp,5);
         }
      }
      //////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////
      // combine trgTPIL, (Xwarp,Ywarp,Zwarp) and subTPIL to obtain (Xout,Yout,Zout)

      Xout = (float *)calloc(trgdim.nv, sizeof(float));
      Yout = (float *)calloc(trgdim.nv, sizeof(float));
      Zout = (float *)calloc(trgdim.nv, sizeof(float));

      combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xwarp, Ywarp, Zwarp, subTPIL);
      combine_warps_and_trans(trgTPIL, Xwarp, Ywarp, Zwarp, Xout, Yout, Zout, trgdim);

      for(int i=0; i<subdim.nv; i++) invXout[i]=invYout[i]=invZout[i]=0.0;

      ivf(Xout, Yout, Zout, trgdim, invXout, invYout, invZout, subdim);

      for(int i=0; i<trgdim.nv; i++) 
      if(trg[i]>0) trg[i]=100;

      tmp=computeReslicedImage2(trg,trgdim,subdim,invXout,invYout,invZout);

      for(int i=0; i<subdim.nv; i++) label[i] += tmp[i]/100.0;

      free(tmp);
      free(trg);
      free(Xout); free(Yout); free(Zout);
   }

   if(opt_v) printf("Thresholding at %6.2f\%\n",threshold);

   for(int i=0; i<subdim.nv; i++)
   {
      label[i] /= (natlas/100.0); // makes label range from 0.0-100.0
      if(label[i] < threshold ) sub[i]=0;
   }

   sprintf(filename,"%s_brainwash.nii",subprefix);
   save_nifti_image(filename, sub, &sub_hdr);

   for(int i=0; i<subdim.nv; i++)
      sub[i] = (short)(label[i]+0.5);
   sprintf(filename,"%s_brainmask.nii",subprefix);
   save_nifti_image(filename, sub, &sub_hdr);

   ////////////////////////////////////////////////////////////////////////////////////////////
   // free allocated memory 
   free(Xwarp); free(Ywarp); free(Zwarp);
   free(invXout); free(invYout); free(invZout);
   free(sub);
   free(trgPIL1); free(trgPIL2); free(trgPIL4); free(trgPIL8);
   free(subPIL1); free(subPIL2); free(subPIL4); free(subPIL8);
   free(mskPIL1); free(mskPIL2); free(mskPIL4); free(mskPIL8);
   ////////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v) printf("THE END\n");
}
