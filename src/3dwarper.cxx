/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  3dwarper.c                                                *
*  Copyright 2004 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <time.h>       
#include <sys/types.h>  
#include <sys/stat.h>  
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include "../include/volume.h"
#include "../include/spm_analyze.h"
#include "../include/babak_lib.h"
#include "../include/sph.h"
#include "../include/smooth.h"
#include "../include/minmax.h"
#include "../include/interpolator.h"

#define YES 1
#define NO 0
#define MAXITER 10

#define MATRIXSIZE 256
#define VOXELSIZE 1.0

extern float detect_lm(SPH &searchsph, SPH &testsph, SHORTIM testim, int lmcm[], SPH &refsph, int lm[]);

void print_matrix(const char * title, float *T);

extern float *resizeXYZ(float *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

extern short *resizeXYZ(short *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

static float v1,v2,v3,v4;
static float w1,w2;

DIM dim1, dim2, dim4, dim8;

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-iter8", 1, '8'},
   {"-iter4", 1, '4'},
   {"-iter2", 1, '2'},
   {"-iter1", 1, '1'},

   {"-trgorient",1,'R'},
   {"-suborient",1,'S'},

   {"-T", 1, 'T'},
   {"-sub", 1, 's'},
   {"-trg", 1, 't'},

   {"-trglm", 1, 'm'},
   {"-sublm", 1, 'M'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h",0,'h'},
   {"-help",0,'h'},
   {"-iter", 1, 'i'},
   {"-u", 1, 'u'},
   {"-o", 1, 'o'},

   {"-thresh", 1, '3'},
   {"-D",0,'D'},
   {"-cubicspline", 0, 'c'},
   {0, 0, 0}
};

int opt_D=NO;
int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit();

short *computeReslicedImage2(short *im1, short *msk, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

short *computeReslicedImage2(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

short *computeReslicedImage2(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

short *computeReslicedImage2(short *,int,int,int,float,float,float,int,int,int,float,float,float,float *,float *,float *,float *);

void print_help_and_exit()
{
   printf("\n\nUsage:\n"
   "\t3dwarper [-v or -verbose] [-h or -help] [-iter N] [-A]\n"
   "\t[-trgorient <orientation code>] [-suborient <orientation code>] [-u <filename>]\n"
   "\t[-o <filename>] [-cubicspline]\n"
   "\t-sub <subject image> -trg <target image>\n\n" 

   "Required arguments:\n\n"
   "\t-sub <subject image>\n"
   "\t\tSpecifies the image to be registred to the <target image>. <subject image> is expected\n"
   "\t\tto be of type short in NIFTI format.\n\n"

   "\t-trg <target image>\n"
   "\t\tSpecifies the image to which the <subject image> is to be registred. <target image> is\n"
   "\t\texpected to be of type short in NIFTI format.\n\n" 

   "Optional arguments:\n\n"
   "\t-v or -verbose\n"
   "\t\tEnables verbose mode.\n\n"

   "\t-h or -help\n"
   "\t\tPrints help message.\n\n"

   "\t-iter N\n"
   "\t\tSpecifies the number of iterations used in finding the initial affine transformation\n"
   "\t\twhen the -A option is specified (default: N=4).\n\n"

   "\t-A\n"
   "\t\tAutomatically finds an initial subject to target affine transformation.\n\n"

   "\t-trgorient <orientation code>\n"
   "\t\tOverrides the orientation information in the target image NIFTI header with the\n"
   "\t\tgiven <orientation code> (PIL, LPS, etc.).\n\n"

   "\t-suborient <orientation code>\n"
   "\t\tOverrides the orientation information in the subject image NIFTI header with the\n"
   "\t\tgiven <orientation code> (PIL, LPS, etc.).\n\n"

   "\t-u <filename>\n"
   "\t\tStores the displacement vector field in the specified <filename>\n"
   "\t\t(default: <filename>=<subject image>_wrp.nii).'\n\n"

   "\t-o <filename>\n"
   "\t\tStores the transformed (registered) <subject image> in the specified <filename>\n"
   "\t\t(default: <filename>=C<subject image>.nii).\n\n"

   "\t-cubicspline\n"
   "\t\tThe output transformed (registered) <subject image> is generated using the\n"
   "\t\tcubic spline interpolation (default is trilinear interpolation).\n\n"
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

void nlReg(short *sub, short *trg, DIM dim, int r, int R, float *Xwarp, float *Ywarp, float *Zwarp, float sd)
{
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
   
   for(int k=0; k<dim.nz; k++)
   for(int j=0; j<dim.ny; j++)
   for(int i=0; i<dim.nx; i++)
   {
      v = k*dim.np + j*dim.nx + i;

      if(trg[v]<=0)
      {
         continue;
      }

      trgsph.set(trgim,i,j,k);

      P[0]=i; P[1]=j; P[2]=k;
      detect_lm(searchsph, subsph, subim, P, trgsph, Q);

      Xw[v] = (Q[0]-P[0])*dim.dx;
      Yw[v] = (Q[1]-P[1])*dim.dy;
      Zw[v] = (Q[2]-P[2])*dim.dz;
   }

   float *Xww, *Yww, *Zww;
   Xww = (float *)calloc(dim.nv, sizeof(float));
   Yww = (float *)calloc(dim.nv, sizeof(float));
   Zww = (float *)calloc(dim.nv, sizeof(float));

   {
      SPH xsph(5);
      SPH ysph(5);
      SPH zsph(5);

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

   tmp=resizeXYZ(Xw,dim.nx,dim.ny,dim.nz,dim.dx,dim.dy,dim.dz,dim1.nx,dim1.ny,dim1.nz,dim1.dx,dim1.dy,dim1.dz);
   free(Xw); Xw=tmp;

   tmp=resizeXYZ(Yw,dim.nx,dim.ny,dim.nz,dim.dx,dim.dy,dim.dz,dim1.nx,dim1.ny,dim1.nz,dim1.dx,dim1.dy,dim1.dz);
   free(Yw); Yw=tmp;

   tmp=resizeXYZ(Zw,dim.nx,dim.ny,dim.nz,dim.dx,dim.dy,dim.dz,dim1.nx,dim1.ny,dim1.nz,dim1.dx,dim1.dy,dim1.dz);
   free(Zw); Zw=tmp;

   for(int i=0; i<dim1.nv; i++) {
      Xwarp[i] += Xw[i];
      Ywarp[i] += Yw[i];
      Zwarp[i] += Zw[i];
   }
   free(Xw); free(Yw); free(Zw);
}

// finds an affine registration in "A" to takes points from sub to trg space
void affineReg(short *sub, short *trg, DIM dim, int r, int R, float *A)
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

      if(sub[v]<=0 || trg[v]<=0)
      {
         flg[v]=0;
         continue;
      }

      subsph.set(subim,i,j,k);

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

   if(opt_v) printMatrix(A,4,4,"",NULL);

   free(PT); 
   free(QT);
   free(flg);
   free(P); 
   free(Q);
}

void generateMultiResolution(short *sub, DIM dim_sub, float *T, float *Xwarp, float *Ywarp, float *Zwarp, short **im1, short **im2, short **im4, short **im8)
{
   float I[16];
   float *Xw,*Yw,*Zw;

   if( *im1 != NULL ) free(*im1);
   if( *im2 != NULL ) free(*im2);
   if( *im4 != NULL ) free(*im4);
   if( *im8 != NULL ) free(*im8);

   Xw = (float *)calloc(dim1.nv, sizeof(float));
   Yw = (float *)calloc(dim1.nv, sizeof(float));
   Zw = (float *)calloc(dim1.nv, sizeof(float));

   for(int i=0; i<dim1.nv; i++) {
      Xw[i] = Xwarp[i];
      Yw[i] = Ywarp[i];
      Zw[i] = Zwarp[i];
   }

   combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xw, Yw, Zw, T);

   *im1=computeReslicedImage2(sub, dim_sub.nx, dim_sub.ny, dim_sub.nz, dim_sub.dx, dim_sub.dy, dim_sub.dz,
   dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xw, Yw, Zw);

   set_to_I(I,4); 
   *im2 = resliceImage(*im1, dim1, dim2, I, LIN); 

   set_to_I(I,4); 
   *im4 = resliceImage(*im2, dim2, dim4, I, LIN); 

   set_to_I(I,4); 
   *im8 = resliceImage(*im4, dim4, dim8, I, LIN); 

   free(Xw); free(Yw); free(Zw);
}

void generateMultiResolution(short *im, DIM im_dim, float *T, short **im1, short **im2, short **im4, short **im8)
{
   float I[16];
   float *invT;		

   if( *im1 != NULL ) free(*im1);
   if( *im2 != NULL ) free(*im2);
   if( *im4 != NULL ) free(*im4);
   if( *im8 != NULL ) free(*im8);

   invT=inv4(T); 
   *im1 = resliceImage(im, im_dim, dim1, invT, LIN); 
   free(invT);

   set_to_I(I,4); 
   *im2 = resliceImage(*im1, dim1, dim2, I, LIN); 

   set_to_I(I,4); 
   *im4 = resliceImage(*im2, dim2, dim4, I, LIN); 

   set_to_I(I,4); 
   *im8 = resliceImage(*im4, dim4, dim8, I, LIN); 
}

void combine_warps_and_trans(float *trgTPIL, float *Xwarp, float *Ywarp, float *Zwarp, float *Xout, float *Yout, float *Zout, DIM dim_trg)
{
   int v;
   float xc,yc,zc;
   float xc1,yc1,zc1;
   float x,y,z;   
   float x1,y1,z1;   
   float i1,j1,k1;   

   xc=dim_trg.dx*(dim_trg.nx-1.0)/2.0;     /* +---+---+ */
   yc=dim_trg.dy*(dim_trg.ny-1.0)/2.0;
   zc=dim_trg.dz*(dim_trg.nz-1.0)/2.0;

   xc1=dim1.dx*(dim1.nx-1.0)/2.0;     /* +---+---+ */
   yc1=dim1.dy*(dim1.ny-1.0)/2.0;
   zc1=dim1.dz*(dim1.nz-1.0)/2.0;
 
   for(int k=0;k<dim_trg.nz;k++) 
   for(int j=0;j<dim_trg.ny;j++) 
   for(int i=0;i<dim_trg.nx;i++) 
   {
      v = k*dim_trg.np + j*dim_trg.nx + i;

      // (i*dx-xc) converts from image coordinates (i,j,z) to (x,y,z) coordinates
      x = (i*dim_trg.dx - xc);
      y = (j*dim_trg.dy - yc);
      z = (k*dim_trg.dz - zc);

      x1 =  trgTPIL[0]*x +trgTPIL[1]*y +trgTPIL[2]*z  +trgTPIL[3];
      y1 =  trgTPIL[4]*x +trgTPIL[5]*y +trgTPIL[6]*z  +trgTPIL[7];
      z1 =  trgTPIL[8]*x +trgTPIL[9]*y +trgTPIL[10]*z +trgTPIL[11];

      i1 = (x1 + xc1) / dim1.dx;
      j1 = (y1 + yc1) / dim1.dy;
      k1 = (z1 + zc1) / dim1.dz;

      x1 += linearInterpolator(i1,j1,k1,Xwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      y1 += linearInterpolator(i1,j1,k1,Ywarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      z1 += linearInterpolator(i1,j1,k1,Zwarp,dim1.nx,dim1.ny,dim1.nz,dim1.np);
      
      Xout[v] = x1 - i*dim_trg.dx + xc;
      Yout[v] = y1 - j*dim_trg.dy + yc;
      Zout[v] = z1 - k*dim_trg.dz + zc;
   }
}

int main(int argc, char **argv)
{
   char affineTransformationFile[1024]="";

   short *subPIL1=NULL, *subPIL2=NULL, *subPIL4=NULL, *subPIL8=NULL;
   short *trgPIL1=NULL, *trgPIL2=NULL, *trgPIL4=NULL, *trgPIL8=NULL;

   FILE *fp;
   char filename[1024];

   char trgprefix[1024]=""; //target image prefix
   char subprefix[1024]=""; //subject image prefix

   getARTHOME();

   nifti_1_header PILhdr;
   sprintf(filename,"%s/PILbrain.nii",ARTHOME);
   PILhdr = read_NIFTI_hdr(filename);

   DIM dim_trg;
   DIM dim_sub;

   dim1.nx = dim1.ny = dim1.nz = MATRIXSIZE;
   dim1.np=dim1.nx*dim1.nx; 
   dim1.nv=dim1.np*dim1.nz; 
   dim1.dx = dim1.dy = dim1.dz = VOXELSIZE;

   dim2.nx = dim2.ny = dim2.nz = MATRIXSIZE/2;
   dim2.np=dim2.nx*dim2.nx; 
   dim2.nv=dim2.np*dim2.nz; 
   dim2.dx = dim2.dy = dim2.dz = VOXELSIZE*2.0;

   dim4.nx = dim4.ny = dim4.nz = MATRIXSIZE/4;
   dim4.np=dim4.nx*dim4.nx; 
   dim4.nv=dim4.np*dim4.nz; 
   dim4.dx = dim4.dy = dim4.dz = VOXELSIZE*4.0;

   dim8.nx = dim8.ny = dim8.nz = MATRIXSIZE/8;
   dim8.np=dim8.nx*dim8.nx; 
   dim8.nv=dim8.np*dim8.nz; 
   dim8.dx = dim8.dy = dim8.dz = VOXELSIZE*8.0;

   // a linear transformation that brings the subject image from its native space to PIL space
   float subTPIL[16]; 

   // a linear transformation that brings the target image from its native space to PIL space
   float trgTPIL[16];

   char trglmfile[1024]="";
   char sublmfile[1024]="";

   int iter8=1;
   int iter4=1;
   int iter2=1;
   int iter1=1;

   nifti_1_header trg_hdr;
   nifti_1_header sub_hdr;
   nifti1_extender extender;

   char subjectImageFile[1024]; 
   char targetImageFile[1024];
   char subOrient[4];  // orientation code for the subject image
   char trgOrient[4];  // orientation code for the target image
   char outputfile[1024]="";
   char warpfile[1024];

   float sub_to_trg[16];
   float *invT;		
   float T[16];	// overall rigid body transformation
   float *Xwarp=NULL, *Ywarp=NULL, *Zwarp=NULL;

   short *sub;
   int Snx,Sny,Snz;
   float Sdx,Sdy,Sdz;

   short *trg;
   int Tnx,Tny,Tnz;
   float Tdx,Tdy,Tdz;

   int thresh=0;

   // number of iterations used in finding the initial affine transformation using -A option
   int niter=4; 

	short *obj;

	int Wx,Wy,Wz;
	int Lx,Ly,Lz;
	int N;	// N=2*L+1

	char *dum;

	int search_win;

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

   /////////////////////////////////////////////////////////////////////
   // important initializations

   warpfile[0]='\0'; 
   subjectImageFile[0]='\0'; 
   targetImageFile[0]='\0'; 
   subOrient[0]='\0'; 
   trgOrient[0]='\0'; 
   /////////////////////////////////////////////////////////////////////

   dum = (char *)malloc(1024);

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
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
         case 'T':
            sprintf(affineTransformationFile,"%s",optarg);
            break;
         case 'R':
            sprintf(trgOrient,"%s",optarg);
            break;
         case 'm':
            sprintf(trglmfile,"%s",optarg);
            break;
         case 'M':
            sprintf(sublmfile,"%s",optarg);
            break;
         case 'S':
            sprintf(subOrient,"%s",optarg);
            break;
         case 'h':
            print_help_and_exit();
         case 'i':
            niter=atoi(optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'D':
            opt_D=YES;
            break;
         case 's':
            sprintf(subjectImageFile,"%s",optarg);
            break;
         case 't':
            sprintf(targetImageFile,"%s",optarg);
            break;
         case '3':
            thresh=atoi(optarg);
            break;
         case 'u':
            sprintf(warpfile,"%s",optarg);
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            break;
         case 'c':
            opt_cubicspline=YES;
            break;
         case '?':
            print_help_and_exit();
		}
	}

   ////////////////////////////////////////////////////////////////////////////////////////////

   // reset unreasonable values to default
   if(iter8 < 0 || iter8>MAXITER ) iter8=1;
   if(iter4 < 0 || iter4>MAXITER ) iter4=1;
   if(iter2 < 0 || iter2>MAXITER ) iter2=1;
   if(iter1 < 0 || iter1>MAXITER ) iter1=1;

   ////////////////////////////////////////////////////////////////////////////////////////////

   // ensure niter is positive and less than MAXITER
   if(niter<0)
   {
      niter=4;
   }
   else if (niter>MAXITER) 
   {
      niter=MAXITER;
   }

   ////////////////////////////////////////////////////////////////////////////////////////////

   // ensure that subject and target images were specified at the command line

   if( subjectImageFile[0] == '\0')
   {
      printf("\nError: you must specify a \"subject image\" using the -sub argument.\n\n");		
      exit(0);
   } 

   if( niftiFilename(subprefix, subjectImageFile)==0 ) { exit(0); }

   if( targetImageFile[0] == '\0')
   {
      printf("\nError: you must specify a \"target image\" using the -trg argument.\n\n");		
      exit(0);
   }

   if( niftiFilename(trgprefix, targetImageFile)==0 ) { exit(0); }

   if(opt_v)
   {
      printf("Subject image file = %s\n",subjectImageFile);
      printf("Target image file = %s\n",targetImageFile);
   }

   sub = (int2 *)read_nifti_image(subjectImageFile, &sub_hdr);
   set_dim(dim_sub,sub_hdr);
   Snx = dim_sub.nx; Sny = dim_sub.ny; Snz = dim_sub.nz;
   Sdx = dim_sub.dx; Sdy = dim_sub.dy; Sdz = dim_sub.dz;

   if(sub==NULL) 
   {
      printf("\nError: Reading subject image %s failed.\n\n",subjectImageFile);
      exit(0);
   }

   trg = (int2 *)read_nifti_image(targetImageFile, &trg_hdr);
   set_dim(dim_trg,trg_hdr);
   Tnx = dim_trg.nx; Tny = dim_trg.ny; Tnz = dim_trg.nz;
   Tdx = dim_trg.dx; Tdy = dim_trg.dy; Tdz = dim_trg.dz;

   if(trg==NULL) 
   {
      printf("\nError: Reading target image %s failed.\n\n",targetImageFile);
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   // find subTPIL 
   if(opt_v) printf("Computing subject image PIL transformation ...\n");
   if(sublmfile[0] != '\0' && opt_v) printf("Subject image landmarks are read from %s\n",sublmfile);
   new_PIL_transform(subjectImageFile, sublmfile, subTPIL);

   // find trgTPIL
   if(opt_v) printf("Computing target image PIL transformation ...\n");
   if(trglmfile[0] != '\0' && opt_v) printf("Target image landmarks are read from %s\n",trglmfile);
   new_PIL_transform(targetImageFile, trglmfile, trgTPIL);

   generateMultiResolution(trg, dim_trg, trgTPIL, &trgPIL1, &trgPIL2, &trgPIL4, &trgPIL8);

   float A[16];

   if( affineTransformationFile[0] != '\0' )
   {
      loadTransformation(affineTransformationFile,sub_to_trg);
      multi(trgTPIL,4,4,sub_to_trg,4,4,subTPIL);  // update subTPIL according to A
   }
   else
   {
      if(opt_v) printf("Affine registration @ 12.5\% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
      affineReg(subPIL8, trgPIL8, dim8, 3, 3, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 25\% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
      affineReg(subPIL4, trgPIL4, dim4, 3, 3, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      if(opt_v) printf("Affine registration @ 50\% resolution ...\n");
      generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
      affineReg(subPIL2, trgPIL2, dim2, 3, 3, A);
      multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

      invT = inv4(trgTPIL);
      multi(invT, 4,4, subTPIL, 4, 4, sub_to_trg);
      free(invT);

      sprintf(filename,"%s_affine.mrx",subprefix);
      fp = fopen(filename,"w");
      printMatrix(sub_to_trg,4,4,"",fp);
      fclose(fp);
   }

   if(opt_v) print_matrix("subject --> target affine transformation",sub_to_trg);

   ////////////////////////////////////////////////////////////////////////////////////////////
   {
      short *tmp;

      invT=inv4(sub_to_trg);
      tmp = resliceImage(sub, dim_sub, dim_trg, invT, LIN); 
      free(invT);

      sprintf(filename,"%s_affine.nii",subprefix);
      if(opt_v) printf("Saving subject image affine transformed to target space in %s ...\n",filename);
      save_nifti_image(filename, tmp, &trg_hdr);
      free(tmp);
   }
   //////////////////////////////////////////////////////////////////
   
   Xwarp = (float *)calloc(dim1.nv, sizeof(float));
   Ywarp = (float *)calloc(dim1.nv, sizeof(float));
   Zwarp = (float *)calloc(dim1.nv, sizeof(float));


   if(opt_v) printf("Non-linear registration @ 12.5\% resolution ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   nlReg(subPIL8, trgPIL8, dim8, 3, 3, Xwarp, Ywarp, Zwarp, 4.0);

   if(opt_v) printf("Non-linear registration @ 25\% resolution ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   nlReg(subPIL4, trgPIL4, dim4, 3, 3, Xwarp, Ywarp, Zwarp, 4.0);

   if(opt_v) printf("Non-linear registration @ 50\% resolution ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   nlReg(subPIL2, trgPIL2, dim2, 3, 3, Xwarp, Ywarp, Zwarp, 4.0);

   if(opt_v) printf("Non-linear registration @ 100\% resolution ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   nlReg(subPIL1, trgPIL1, dim1, 3, 3, Xwarp, Ywarp, Zwarp, 2.0);


   generateMultiResolution(sub, dim_sub, subTPIL, Xwarp, Ywarp, Zwarp, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   set_dim(PILhdr,dim1);
   sprintf(filename,"%s_PIL1.nii",trgprefix);
   save_nifti_image(filename, trgPIL1, &PILhdr);
   sprintf(filename,"%s_PIL1.nii",subprefix);
   save_nifti_image(filename, subPIL1, &PILhdr);

   {
      short *Csub;
      float *Xout, *Yout, *Zout;
      float T[16];

      Xout = (float *)calloc(dim_trg.nv, sizeof(float));
      Yout = (float *)calloc(dim_trg.nv, sizeof(float));
      Zout = (float *)calloc(dim_trg.nv, sizeof(float));

      combine_warps_and_trans(dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim1.dy, dim1.dz, Xwarp, Ywarp, Zwarp, subTPIL);

      combine_warps_and_trans(trgTPIL, Xwarp, Ywarp, Zwarp, Xout, Yout, Zout, dim_trg);

      Csub=computeReslicedImage2(sub, dim_sub, dim_trg, Xout, Yout, Zout);

      if(outputfile[0]=='\0') sprintf(filename,"C%s.nii",subprefix);
      else sprintf(filename,"%s",outputfile);
      if(opt_v) printf("Saving output (warped) image file %s ...\n",filename);
      save_nifti_image(filename, Csub, &trg_hdr);

      free(Xout); free(Yout); free(Zout);
      free(Csub);
   }

exit(0);
   ////////////////////////////////////////////////////////////////////////////////////////////

   {
      float min, max, s=0.0;
      short *sdum;

      //////////////////////////////////////////////////////////////////////////////////////////
      /// If no displacement field filename is specified, the program automatically sets it to be
      /// the subject filename appended by '_wrp'. For example, if the subject filename is
      /// test.nii, the warp parameters filename will be test_wrp.nii
      //////////////////////////////////////////////////////////////////////////////////////////
      if(warpfile[0]=='\0') 
      {
         // returns the filename in dum, without the path information.
         getfilename(dum,subjectImageFile);

         sprintf(warpfile,"%s_wrp.nii", strsep(&dum, "."));
      }

      if(opt_v)
      {
         printf("\nOutput warp parameters file = %s\n",warpfile);
      }
      //////////////////////////////////////////////////////////////////////////////////////////

      combine_warps_and_trans(Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, sub_to_trg);
 
      minmax(Xwarp, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Ywarp, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Zwarp, dim_trg.nv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      fp=fopen(warpfile,"w");
      trg_hdr.bitpix=8*sizeof(short);
      trg_hdr.vox_offset=352.0;
      trg_hdr.dim[0]=5;
      trg_hdr.dim[4]=1;
      trg_hdr.dim[5]=3;
      trg_hdr.pixdim[4] = 0.0;
      trg_hdr.datatype=DT_SIGNED_SHORT;
      trg_hdr.intent_code=NIFTI_INTENT_VECTOR;
      trg_hdr.scl_slope=s/32000;
      extender.extension[0]=0;
      fwrite(&trg_hdr,sizeof(nifti_1_header),1,fp);
      fwrite(&extender, sizeof(nifti1_extender),1,fp);

      sdum = (short *)calloc(dim_trg.nv, sizeof(short));
      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Xwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Ywarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      for(int i=0; i<dim_trg.nv; i++)
      {
         sdum[i] = (short)(Zwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),dim_trg.nv,fp);

      fclose(fp);

      free(sdum);
   }

   free(Xwarp); free(Ywarp); free(Zwarp);
   free(sub);
   free(trg);

   if(opt_v) printf("THE END\n");
}

short *computeReslicedImage2(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp)
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
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

	np1=nx1*ny1;

	im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	for(int j=0;j<ny2;j++) 
  	for(int i=0;i<nx2;i++) 
	{
		z = (k*dz2 - zc2 + Zwarp[q] + zc1) /dz1;
		y = (j*dy2 - yc2 + Ywarp[q] + yc1) /dy1;
		x = (i*dx2 - xc2 + Xwarp[q] + xc1) /dx1;

		if(opt_cubicspline)
			im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
		else
			im2[q++]=(short)(linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

short *computeReslicedImage2(short *im1, short *msk, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp)
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
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

	np1=nx1*ny1;

	im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	for(int j=0;j<ny2;j++) 
  	for(int i=0;i<nx2;i++) 
	{
		if(!msk[q])
		{
			im2[q++]=0;
		}
		else
		{
			z = (k*dz2 - zc2 + Zwarp[q] + zc1) /dz1;
			y = (j*dy2 - yc2 + Ywarp[q] + yc1) /dy1;
       		x = (i*dx2 - xc2 + Xwarp[q] + xc1) /dx1;

			if(opt_cubicspline)
				im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
			else
				im2[q++] = (short)(linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
		}
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

short *computeReslicedImage2(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp)
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
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

	np1=nx1*ny1;

   im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));
   if(im2==NULL)
   {
      printf("Error: memory allocation failure for variable 'im2', aborting ...\n");
      exit(1);
   }

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	for(int j=0;j<ny2;j++) 
  	for(int i=0;i<nx2;i++) 
	{

		z = (k*dz2 - zc2 + Zwarp[q] + zc1) /dz1;
		y = (j*dy2 - yc2 + Ywarp[q] + yc1) /dy1;
		x = (i*dx2 - xc2 + Xwarp[q] + xc1) /dx1;

			if(opt_cubicspline)
            {
				im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
            }
			else
            {
	   			im2[q++] = (short)(linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
            }
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

short *computeReslicedImage2(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
  	float  x,y,z;   
  	float  xx,yy,zz;   
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
		c = (float *)calloc(nx1*ny1*nz1, sizeof(float));
		cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
	}

	np1=nx1*ny1;

	im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
	for(int k=0;k<nz2;k++) 
	{
		for(int j=0;j<ny2;j++) 
		{
  			for(int i=0;i<nx2;i++) 
			{
				zz = k*dz2 - zc2 + Zwarp[q];
				yy = j*dy2 - yc2 + Ywarp[q];
				xx = i*dx2 - xc2 + Xwarp[q];

				x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3]   + xc1 )/dx1;
				y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7]   + yc1 )/dy1;
				z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]  + zc1 )/dz1;

				if(opt_cubicspline)
					im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
				else
					im2[q++] = (short)(linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
			}
		}
	}

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

void print_matrix(const char * title, float *T)
{
	printf("\n%s:",title);
	printf("\n%f\t%f\t%f\t%f",T[0],T[1],T[2],T[3]);
	printf("\n%f\t%f\t%f\t%f",T[4],T[5],T[6],T[7]);
	printf("\n%f\t%f\t%f\t%f",T[8],T[9],T[10],T[11]);
	printf("\n%f\t%f\t%f\t%f",T[12],T[13],T[14],T[15]);
	printf("\n");
}

