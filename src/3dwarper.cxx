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

   {"-sub", 1, 's'},
   {"-trg", 1, 't'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},
   {"-h",0,'h'},
   {"-help",0,'h'},
   {"-iter", 1, 'i'},
   {"-u", 1, 'u'},
   {"-o", 1, 'o'},

   {"-thresh", 1, '3'},
   {"-w", 1, '7'},
   {"-s", 1, '9'},
   {"-D",0,'D'},
   {"-sd", 1, 'd'},
   {"-cubicspline", 0, 'c'},
   {0, 0, 0}
};

int opt_D=NO;
int opt_w=NO;
int opt_s=NO;
int opt_sd=NO;
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
   "\t[-o <filename>] [-cubicspline] [-sd N] [-w N] [-s N]\n"
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

   "\t-sd N\n"
   "\t\tSpecifies the degree of smoothing applied to the displacement vector field.\n"
   "\t\tN typically ranges from 5.0 to 12.0 mm (default = search window size minus 1 mm).\n\n"

   "\t-w N\n"
   "\t\tCorrelation window size in voxels (default=5; minimum=3).\n\n"

   "\t-s N\n"
   "\t\tSearch window size in voxels (default=5; minimum=3).\n\n"
	);
	exit(0);
}

void computeWarpField(
float resFactor, 
short *obj, int Onx, int Ony, int Onz,float Odx, float Ody, float Odz, 
short *trg, int Tnx, int Tny, int Tnz,float Tdx, float Tdy, float Tdz, 
float *Xwarp, float *Ywarp, float *Zwarp,
int Lx, int Wx,
float sd, int thresh)
{
	int Ly, Lz;
	int Wy, Wz;
	int Tnv;
	float V;
	float *ARobj;	// array extracted from object and target images
	float *ARtrg;

	short *HRtrg;	// high res. target and object images
	short *HRobj;
	int HRnx, HRny, HRnz, HRnp;
	float HRdx, HRdy, HRdz;
	float *Xw, *Yw, *Zw;

	int xopt, yopt, zopt;
	int N3;	// N*N*N
	float CC, CCMAX;
	float Sx, Sx2, Sy, Sxy;
	float num,den;

	float *Xwarp_tmp, *Ywarp_tmp, *Zwarp_tmp;

	Ly = (int)(Lx*Tdx/Tdy + 0.5); if(Ly==0) Ly++;
	Lz = (int)(Lx*Tdx/Tdz + 0.5); if(Lz==0) Lz++;
	N3 = (2*Lx+1)*(2*Ly+1)*(2*Lz+1);

	Wy = (int)(Wx*Tdx/Tdy + 0.5); if(Wy==0) Wy++;
	Wz = (int)(Wx*Tdx/Tdz + 0.5); if(Wz==0) Wz++;

	Tnv=Tnx*Tny*Tnz;

	ARobj=(float *)calloc(N3,sizeof(float));
	ARtrg=(float *)calloc(N3,sizeof(float));

	HRdx = (float)(	Tdx * resFactor );
	HRdy = (float)( Tdy * resFactor );
	HRdz = (float)( Tdz * resFactor );

	HRnx = (int)(Tnx/resFactor + 0.5 );
	HRny = (int)(Tny/resFactor + 0.5 );
	HRnz = (int)(Tnz/resFactor + 0.5 );

   if(opt_v)
   {
      printf("\nResolution reduction factor = %4.2f",resFactor);
      printf("\nMatrix size = %d x %d x %d (voxels)", HRnx, HRny, HRnz);
      printf("\nVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", HRdx,HRdy,HRdz);
   }

	HRnp = HRnx * HRny;

	// At first iteration, Xwapr, Ywarp, Zwarp are zero
	Xw=resizeXYZ(Xwarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
	Yw=resizeXYZ(Ywarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
	Zw=resizeXYZ(Zwarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);

   {
      float *tmp=NULL;
      float sdx, sdy, sdz;

      // added 5/11/11 due to bug report by Tito
      // have to make sure that x>0.0 in sqrt(x)
      if( (HRdx*HRdx - Odx*Odx) < 0.0 )
      {
         sdx=0.0;
      }
      else
      {
         sdx=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - Odx*Odx ) )/Odx );
      }

      if( (HRdy*HRdy - Ody*Ody) < 0.0 )
      {
         sdy=0.0;
      }
      else
      {
         sdy=(float)( sqrt( (0.5/log(2.0)) * ( HRdy*HRdy - Ody*Ody ) )/Ody );
      }

      if( (HRdz*HRdz - Odz*Odz) < 0.0 )
      {
         sdz=0.0;
      }
      else
      {
         sdz=(float)( sqrt( (0.5/log(2.0)) * ( HRdz*HRdz - Odz*Odz ) )/Odz );
      }

      tmp = smoothXYZ(obj,Onx,Ony,Onz,sdx,sdy,sdz);

      if(tmp==NULL)
      {
         printf("Error: smoothXYZ() failed to return an image! Aborting ...\n");
         exit(1);
      }

      HRobj=computeReslicedImage2(tmp, Onx, Ony, Onz, Odx, Ody, Odz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Xw, Yw, Zw);

      free(tmp);
   }

	HRtrg=resizeXYZ(trg, Tnx ,Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);

	for(int k=0; k<HRnz; k++)
	for(int j=0; j<HRny; j++)
	for(int i=0; i<HRnx; i++)
	{
		if(HRtrg[k*HRnp + j*HRnx + i]<=thresh)
		{
			Xw[k*HRnp + j*HRnx + i] = 0.0;
			Yw[k*HRnp + j*HRnx + i] = 0.0;
			Zw[k*HRnp + j*HRnx + i] = 0.0;
			continue;
		}

		extractArray(HRtrg, HRnx, HRny, HRnz, i, j, k, Lx, Ly, Lz, ARtrg);

		Sy=0.0;
		for(int n=0; n<N3; n++) Sy += ARtrg[n];

		if( Sy == 0.0 )
		{
			Xw[k*HRnp + j*HRnx + i] = 0.0;
			Yw[k*HRnp + j*HRnx + i] = 0.0;
			Zw[k*HRnp + j*HRnx + i] = 0.0;
			continue;
		}

		CCMAX=0.0; 	// IMPORTANT: we are not interested in -tive correlations
					// if CMAX is set to -1, program gives unexpected results

		xopt=yopt=zopt=0;

		for(int x=-Wx; x<=Wx; x++)
		for(int y=-Wy; y<=Wy; y++)
		for(int z=-Wz; z<=Wz; z++)
		{
			extractArray(HRobj, HRnx, HRny, HRnz, i+x, j+y, k+z, Lx, Ly, Lz, ARobj);

			Sx=Sx2=Sxy=0.0;
			for(int n=0; n<N3; n++)
			{
				V = ARobj[n];
				Sx += V;
				Sx2 += (V*V);
				Sxy += (V*ARtrg[n]);
			}

			num = Sxy-Sx*Sy/N3;
			den = (float)sqrt( (double)(Sx2-Sx*Sx/N3) );

			if(den==0.0) continue;

			CC = num/den;

			if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; zopt=z; }
		}

		Xw[k*HRnp + j*HRnx + i] = xopt*HRdx;
		Yw[k*HRnp + j*HRnx + i] = yopt*HRdy;
		Zw[k*HRnp + j*HRnx + i] = zopt*HRdz;
	}

	medianFilter(Xw, HRnx, HRny, HRnz, Wx, Wy, Wz);
	medianFilter(Yw, HRnx, HRny, HRnz, Wx, Wy, Wz);
	medianFilter(Zw, HRnx, HRny, HRnz, Wx, Wy, Wz);

	Xwarp_tmp=resizeXYZ(Xw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
	Ywarp_tmp=resizeXYZ(Yw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
	Zwarp_tmp=resizeXYZ(Zw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
	free(Xw); free(Yw); free(Zw);

	for(int n=0; n<Tnv; n++)
	{
		Xwarp[n] += Xwarp_tmp[n];
		Ywarp[n] += Ywarp_tmp[n];
		Zwarp[n] += Zwarp_tmp[n];
	}
	free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

	if(resFactor==1.0) sd = Wx;
  
	Xwarp_tmp=smoothXYZ(Xwarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);
	Ywarp_tmp=smoothXYZ(Ywarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);
	Zwarp_tmp=smoothXYZ(Zwarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);

	for(int i=0; i<Tnv; i++)
	{
		Xwarp[i] = Xwarp_tmp[i];
		Ywarp[i] = Ywarp_tmp[i];
		Zwarp[i] = Zwarp_tmp[i];
	}

	free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

	free(HRobj);
	free(HRtrg);

	free(ARobj);
	free(ARtrg);
}

void computeWarpField(float resFactor, short *obj, DIM O, short *trg, DIM T, float *Xwarp, float *Ywarp, float *Zwarp,
int Lx, int Wx, float sd, int thresh)
{
	int Ly, Lz;
	int Wy, Wz;
	int Tnv;
	float V;
	float *ARobj;	// array extracted from object and target images
	float *ARtrg;

	short *HRtrg;	// high res. target and object images
	short *HRobj;
	int HRnx, HRny, HRnz, HRnp;
	float HRdx, HRdy, HRdz;
	float *Xw, *Yw, *Zw;

	int xopt, yopt, zopt;
	int N3;	// N*N*N
	float CC, CCMAX;
	float Sx, Sx2, Sy, Sxy;
	float num,den;

	float *Xwarp_tmp, *Ywarp_tmp, *Zwarp_tmp;

	Ly = (int)(Lx*T.dx/T.dy + 0.5); if(Ly==0) Ly++;
	Lz = (int)(Lx*T.dx/T.dz + 0.5); if(Lz==0) Lz++;
	N3 = (2*Lx+1)*(2*Ly+1)*(2*Lz+1);

	Wy = (int)(Wx*T.dx/T.dy + 0.5); if(Wy==0) Wy++;
	Wz = (int)(Wx*T.dx/T.dz + 0.5); if(Wz==0) Wz++;

	Tnv=T.nx*T.ny*T.nz;

	ARobj=(float *)calloc(N3,sizeof(float));
	ARtrg=(float *)calloc(N3,sizeof(float));

	HRdx = (float)(	T.dx * resFactor );
	HRdy = (float)( T.dy * resFactor );
	HRdz = (float)( T.dz * resFactor );

	HRnx = (int)(T.nx/resFactor + 0.5 );
	HRny = (int)(T.ny/resFactor + 0.5 );
	HRnz = (int)(T.nz/resFactor + 0.5 );

   if(opt_v)
   {
      printf("\nResolution reduction factor = %4.2f",resFactor);
      printf("\nMatrix size = %d x %d x %d (voxels)", HRnx, HRny, HRnz);
      printf("\nVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", HRdx,HRdy,HRdz);
   }

	HRnp = HRnx * HRny;

	// At first iteration, Xwapr, Ywarp, Zwarp are zero
	Xw=resizeXYZ(Xwarp, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
	Yw=resizeXYZ(Ywarp, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
	Zw=resizeXYZ(Zwarp, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);

   {
      float *tmp;
      float sdx, sdy, sdz;

      // edited added 5/11/11 due to bug report by Tito
      // have to make sure that x>0.0 in sqrt(x)
      if( (HRdx*HRdx - O.dx*O.dx) < 0.0 )
      {
         sdx=0.0;
      }
      else
      {
         sdx=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - O.dx*O.dx ) )/O.dx );
      }

      if( (HRdy*HRdy - O.dy*O.dy) < 0.0 )
      {
         sdy=0.0;
      }
      else
      {
         sdy=(float)( sqrt( (0.5/log(2.0)) * ( HRdy*HRdy - O.dy*O.dy ) )/O.dy );
      }

      if( (HRdz*HRdz - O.dz*O.dz) < 0.0 )
      {
         sdz=0.0;
      }
      else
      {
         sdz=(float)( sqrt( (0.5/log(2.0)) * ( HRdz*HRdz - O.dz*O.dz ) )/O.dz );
      }

      tmp = smoothXYZ(obj,O.nx,O.ny,O.nz,sdx,sdy,sdz);

      if(tmp==NULL)
      {
         printf("Error: smoothXYZ() failed to return an image! Aborting ...\n");
         exit(1);
      }

      HRobj=computeReslicedImage2(tmp, O.nx, O.ny, O.nz, O.dx, O.dy, O.dz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Xw, Yw, Zw);
      free(tmp);
   }

	HRtrg=resizeXYZ(trg, T.nx ,T.ny, T.nz, T.dx, T.dy, T.dz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);

	for(int k=0; k<HRnz; k++)
	for(int j=0; j<HRny; j++)
	for(int i=0; i<HRnx; i++)
	{
		if(HRtrg[k*HRnp + j*HRnx + i]<=thresh)
		{
			Xw[k*HRnp + j*HRnx + i] = 0.0;
			Yw[k*HRnp + j*HRnx + i] = 0.0;
			Zw[k*HRnp + j*HRnx + i] = 0.0;
			continue;
		}

		extractArray(HRtrg, HRnx, HRny, HRnz, i, j, k, Lx, Ly, Lz, ARtrg);

		Sy=0.0;
		for(int n=0; n<N3; n++) Sy += ARtrg[n];

		if( Sy == 0.0 )
		{
			Xw[k*HRnp + j*HRnx + i] = 0.0;
			Yw[k*HRnp + j*HRnx + i] = 0.0;
			Zw[k*HRnp + j*HRnx + i] = 0.0;
			continue;
		}

		CCMAX=0.0; 	// IMPORTANT: we are not interested in -tive correlations
					// if CMAX is set to -1, program given unexpected results

		xopt=yopt=zopt=0;

		for(int x=-Wx; x<=Wx; x++)
		for(int y=-Wy; y<=Wy; y++)
		for(int z=-Wz; z<=Wz; z++)
		{
			extractArray(HRobj, HRnx, HRny, HRnz, i+x, j+y, k+z, Lx, Ly, Lz, ARobj);

			Sx=Sx2=Sxy=0.0;
			for(int n=0; n<N3; n++)
			{
				V = ARobj[n];
				Sx += V;
				Sx2 += (V*V);
				Sxy += (V*ARtrg[n]);
			}

			num = Sxy-Sx*Sy/N3;
			den = (float)sqrt( (double)(Sx2-Sx*Sx/N3) );

			if(den==0.0) continue;

			CC = num/den;

			if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; zopt=z; }
		}

		Xw[k*HRnp + j*HRnx + i] = xopt*HRdx;
		Yw[k*HRnp + j*HRnx + i] = yopt*HRdy;
		Zw[k*HRnp + j*HRnx + i] = zopt*HRdz;
	}

	medianFilter(Xw, HRnx, HRny, HRnz, Wx, Wy, Wz);
	medianFilter(Yw, HRnx, HRny, HRnz, Wx, Wy, Wz);
	medianFilter(Zw, HRnx, HRny, HRnz, Wx, Wy, Wz);

	Xwarp_tmp=resizeXYZ(Xw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz);
	Ywarp_tmp=resizeXYZ(Yw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz);
	Zwarp_tmp=resizeXYZ(Zw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, T.nx, T.ny, T.nz, T.dx, T.dy, T.dz);
	free(Xw); free(Yw); free(Zw);

	for(int n=0; n<Tnv; n++)
	{
		Xwarp[n] += Xwarp_tmp[n];
		Ywarp[n] += Ywarp_tmp[n];
		Zwarp[n] += Zwarp_tmp[n];
	}
	free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

	if(resFactor==1.0) sd = Wx;

	Xwarp_tmp=smoothXYZ(Xwarp, T.nx, T.ny, T.nz, sd, sd*Wy/Wx, sd*Wz/Wx);
	Ywarp_tmp=smoothXYZ(Ywarp, T.nx, T.ny, T.nz, sd, sd*Wy/Wx, sd*Wz/Wx);
	Zwarp_tmp=smoothXYZ(Zwarp, T.nx, T.ny, T.nz, sd, sd*Wy/Wx, sd*Wz/Wx);

	for(int i=0; i<Tnv; i++)
	{
		Xwarp[i] = Xwarp_tmp[i];
		Ywarp[i] = Ywarp_tmp[i];
		Zwarp[i] = Zwarp_tmp[i];
	}

	free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

	free(HRobj);
	free(HRtrg);

	free(ARobj);
	free(ARtrg);
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

   if(opt_v) printMatrix(A,4,4,"A:",NULL);

   free(PT); 
   free(QT);
   free(flg);
   free(P); 
   free(Q);
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

int main(int argc, char **argv)
{
   FILE *fp;
   char filename[1024];

   char trgprefix[1024]=""; //target image prefix
   char subprefix[1024]=""; //subject image prefix

   getARTHOME();

   nifti_1_header PILhdr;
   sprintf(filename,"%s/PILbrain.nii",ARTHOME);
   PILhdr = read_NIFTI_hdr(filename);

   short *subPIL1=NULL, *subPIL2=NULL, *subPIL4=NULL, *subPIL8=NULL;
   short *trgPIL1=NULL, *trgPIL2=NULL, *trgPIL4=NULL, *trgPIL8=NULL;

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

   char sublmfile[1024]="";
   char trglmfile[1024]="";

   int iter8=1;
   int iter4=1;
   int iter2=1;
   int iter1=1;

   nifti_1_header trg_hdr;
   nifti1_extender extender;

   char subjectImageFile[1024]; 
   char targetImageFile[1024];
   char subOrient[4];  // orientation code for the subject image
   char trgOrient[4];  // orientation code for the target image
   char outputfile[1024];
   char warpfile[1024];

   float sub_to_trg[16];
   float *invT;		
   float T[16];		// overall rigid body transformation
   float *Xwarp, *Ywarp, *Zwarp;

   short *sub;
   int Snx,Sny,Snz,Snv;
   float Sdx,Sdy,Sdz;

   short *trg;
   int Tnx,Tny,Tnz,Tnv;
   float Tdx,Tdy,Tdz;

   int thresh=0;

   // number of iterations used in finding the initial affine transformation using -A option
   int niter=4; 

	short *obj;

	float sd;
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

   outputfile[0]='\0'; 
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
         case 'R':
            sprintf(trgOrient,"%s",optarg);
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
         case '7':
            N=atoi(optarg);
            opt_w=YES;
            break;
         case '9':
            search_win=atoi(optarg);
            opt_s=YES;
            break;
         case 'd':
            sd = atof(optarg);
            opt_sd=YES;
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

   sub=readNiftiImage(subjectImageFile, &dim_sub, opt_v);
   Snv=dim_sub.nx*dim_sub.ny*dim_sub.nz;
   Snx = dim_sub.nx; Sny = dim_sub.ny; Snz = dim_sub.nz;
   Sdx = dim_sub.dx; Sdy = dim_sub.dy; Sdz = dim_sub.dz;

   if(sub==NULL) 
   {
      printf("\nError: Reading subject image %s failed.\n\n",subjectImageFile);
      exit(0);
   }

   trg=readNiftiImage(targetImageFile, &dim_trg, opt_v);
   Tnx = dim_trg.nx; Tny = dim_trg.ny; Tnz = dim_trg.nz;
   Tdx = dim_trg.dx; Tdy = dim_trg.dy; Tdz = dim_trg.dz;
   Tnv=Tnx*Tny*Tnz;

   if(trg==NULL) 
   {
      printf("\nError: Reading target image %s failed.\n\n",targetImageFile);
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   // find subTPIL 
   if(opt_v) printf("Computing subject image PIL transformation ...\n");
   new_PIL_transform(subjectImageFile, sublmfile, subTPIL);

   // find trgTPIL
   if(opt_v) printf("Computing target image PIL transformation ...\n");
   new_PIL_transform(targetImageFile, trglmfile, trgTPIL);

   generateMultiResolution(trg, dim_trg, trgTPIL, &trgPIL1, &trgPIL2, &trgPIL4, &trgPIL8);

   float A[16];
   if(opt_v) printf("Processing resolution 1/8 ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   affineReg(subPIL8, trgPIL8, dim8, 3, 3, A);
   multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

   if(opt_v) printf("Processing resolution 1/4 ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   affineReg(subPIL4, trgPIL4, dim4, 3, 3, A);
   multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

   if(opt_v) printf("Processing resolution 1/2 ...\n");
   generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);
   affineReg(subPIL2, trgPIL2, dim2, 3, 3, A);
   multi(A,4,4,subTPIL,4,4,subTPIL);  // update subTPIL according to A

   invT = inv4(trgTPIL);
   multi(invT, 4,4, subTPIL, 4, 4, sub_to_trg);
   free(invT);

   if(opt_v) print_matrix("subject --> target affine transformation",sub_to_trg);

   sprintf(filename,"%s_affine.mrx",subprefix);
   fp = fopen(filename,"w");
   printMatrix(sub_to_trg,4,4,"",fp);
   fclose(fp);

   ////////////////////////////////////////////////////////////////////////////////////////////
   {
      short *tmp;

      trg_hdr = read_NIFTI_hdr(targetImageFile);
      invT=inv4(sub_to_trg);
      tmp = resliceImage(sub, dim_sub, dim_trg, invT, LIN); 
      free(invT);

      sprintf(filename,"%s_affine.nii",subprefix);
      save_nifti_image(filename, tmp, &trg_hdr);
      free(tmp);
   }
   //////////////////////////////////////////////////////////////////
/*
   generateMultiResolution(sub, dim_sub, subTPIL, &subPIL1, &subPIL2, &subPIL4, &subPIL8);

   set_dim(PILhdr,dim1);
   sprintf(filename,"%s_PIL1.nii",trgprefix);
   save_nifti_image(filename, trgPIL1, &PILhdr);
   sprintf(filename,"%s_PIL1.nii",subprefix);
   save_nifti_image(filename, subPIL1, &PILhdr);

   set_dim(PILhdr,dim2);
   sprintf(filename,"%s_PIL2.nii",trgprefix);
   save_nifti_image(filename, trgPIL2, &PILhdr);
   sprintf(filename,"%s_PIL2.nii",subprefix);
   save_nifti_image(filename, subPIL2, &PILhdr);

   set_dim(PILhdr,dim4);
   sprintf(filename,"%s_PIL4.nii",trgprefix);
   save_nifti_image(filename, trgPIL4, &PILhdr);
   sprintf(filename,"%s_PIL4.nii",subprefix);
   save_nifti_image(filename, subPIL4, &PILhdr);

   set_dim(PILhdr,dim8);
   sprintf(filename,"%s_PIL8.nii",trgprefix);
   save_nifti_image(filename, trgPIL8, &PILhdr);
   sprintf(filename,"%s_PIL8.nii",subprefix);
   save_nifti_image(filename, subPIL8, &PILhdr);
*/
   ////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////

   if(!opt_w || N<=0) N=5;
   if( (N%2)==0 ) N+=1;			// make sure it's odd
   if(N==1) N=3;

   if(!opt_s || search_win<=0) search_win=5;
   if( (search_win%2)==0 ) search_win+=1;			// make sure it's odd
   if(search_win==1) search_win=3;

   Lx = (N-1)/2; if(Lx==0) Lx++;
   Ly = (int)(Lx*Tdx/Tdy + 0.5); if(Ly==0) Ly++;
   Lz = (int)(Lx*Tdx/Tdz + 0.5); if(Lz==0) Lz++;

   Wx = (search_win-1)/2; if(Wx==0) Wx++;
   Wy = (int)(Wx*Tdx/Tdy + 0.5); if(Wy==0) Wy++;
   Wz = (int)(Wx*Tdx/Tdz + 0.5); if(Wz==0) Wz++;

   if(!opt_sd || sd<0.0)
      sd = 2.0*Wx;

   ////////////////////////////////////////////////////////////////////////////////////////////

   // allocate memory for and initialize warp field
   Xwarp=(float *)calloc(Tnv,sizeof(float));
   Ywarp=(float *)calloc(Tnv,sizeof(float));
   Zwarp=(float *)calloc(Tnv,sizeof(float));

   if(Xwarp==NULL || Ywarp==NULL || Zwarp==NULL) 
   {
      printf("\n\nMemory allocation error, aborting ...\n\n");
      exit(0);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////


   //////////////////////////////////////////////////////////////////
   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
      printf("Non-linear registration ...\n");
   }

   for(int n=0; n<Tnv; n++) Xwarp[n]=Ywarp[n]=Zwarp[n]=0.0;

   for(int i=0; i<iter8; i++)
   {
      computeWarpField(8.0, obj, Snx, Sny, Snz, Sdx, Sdy, Sdz, trg, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, Lx, Wx, sd, thresh);
   }

   for(int i=0; i<iter4; i++)
   {
      computeWarpField(4.0, obj, Snx, Sny, Snz, Sdx, Sdy, Sdz, trg, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, Lx, Wx, sd, thresh);
   }

   for(int i=0; i<iter2; i++)
   {
      computeWarpField(2.0, obj, Snx, Sny, Snz, Sdx, Sdy, Sdz, trg, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, Lx, Wx, sd, thresh);
   }

   for(int i=0; i<iter1; i++)
   {
      computeWarpField(1.0, obj, Snx, Sny, Snz, Sdx, Sdy, Sdz, trg, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, Lx, Wx, sd, thresh);
   }

   free(obj);

   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
   }
   //////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////
   // save warped subject image

   invT=inv4(sub_to_trg);
   obj=computeReslicedImage2(sub, Snx, Sny, Snz, Sdx, Sdy, Sdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, Xwarp, Ywarp, Zwarp, invT);
   free(invT);
		
   //////////////////////////////////////////////////////////////////////////////////////////
   /// If no output filename is specified, the program automatically sets it by putting the
   /// letter 'C' in front of the subject filename. For example, if the subject filename is 
   /// test.nii, the output filename will be Ctest.nii.
   //////////////////////////////////////////////////////////////////////////////////////////
   if(outputfile[0]=='\0') 
   {
      // returns the filename in dum, without the path information.
      getfilename(dum,subjectImageFile); 

      sprintf(outputfile,"C%s",dum);
   }

   if(opt_v)
   {
      printf("\nOutput (warped) image file = %s\n",outputfile);
   }

   save_nifti_image(outputfile, obj, &trg_hdr);

   free(obj);
   //////////////////////////////////////////////////////////////////

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
 
      minmax(Xwarp, Tnv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Ywarp, Tnv, min, max);
      if(max>s) s=max;
      if(-min>s) s=-min;

      minmax(Zwarp, Tnv, min, max);
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

      sdum = (short *)calloc(Tnv, sizeof(short));
      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Xwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Ywarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

      for(int i=0; i<Tnv; i++)
      {
         sdum[i] = (short)(Zwarp[i]*32000/s + 0.5);
      }
      fwrite(sdum,sizeof(short),Tnv,fp);

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

