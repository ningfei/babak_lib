#define _smooth

#include <stdlib.h>
#include <stdio.h>
#include "babak_lib.h"

// NOTE: The standard deviation (SD) in all these functions is in units of pixels

/////////////////////////////////////////////////////////////////////
float *smoothX(short *image, int nx, int ny, float sd)
{
   int n;
   int np;
   float *h;
   float *image2;
   short *x;

   np=nx*ny;

   image2=(float *)calloc(np,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd,&n);

   for(int j=0;j<ny;j++)
   {
      x = image + nx*j;

      for(int i=0;i<nx;i++) 
      {
         image2[nx*j+i]=conv_pnt_sk(x,nx,h,n,i);
      }
   }

   free(h);

   return(image2);
}

float *smoothX(float *image, int nx, int ny, float sd)
{
   int n;
   int np;
   float *h;
   float *image2;
   float *x;

   np=nx*ny;

   image2=(float *)calloc(np,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd,&n);

   for(int j=0;j<ny;j++)
   {
      for(int i=0;i<nx;i++)
      {
         x=image + nx*j;

         image2[nx*j+i]=conv_pnt_sk(x,nx,h,n,i);
      }
   }

   free(h);

   return(image2);
}

float *smoothX(short *image, int nx, int ny, int nz, float sd)
{
   int n;
   int slice_offset;
   int row_offset;

   float *h;
   float *image2;
   short *x;

   image2=(float *)calloc(nx*ny*nz,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   for(int k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int j=0;j<ny;j++)
      {
         row_offset = slice_offset + nx*j;
         for(int i=0;i<nx;i++)
         {
            x = image + row_offset;
            image2[ row_offset + i ] = conv_pnt_sk(x,nx,h,n,i);
         }
      }
   }

   free(h);

   return(image2);
}

float *smoothX(float *image, int nx, int ny, int nz, float sd)
{
   int n;
   int slice_offset;
   int row_offset;

   float *h;
   float *image2;
   float *x;

   image2=(float *)calloc(nx*ny*nz,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   for(int k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int j=0;j<ny;j++)
      {
         row_offset = slice_offset + nx*j;
         for(int i=0;i<nx;i++)
         {
            x = image + row_offset;
            image2[ row_offset + i ] = conv_pnt_sk(x,nx,h,n,i);
         }
      }
   }

   free(h);

   return(image2);
}

/////////////////////////////////////////////////////////////////////

float *smoothY(float *image, int nx, int ny, float sd)
{
   int n;
   int np;
   float *h;
   float *image2;
   float *y;

   np=nx*ny;

   image2=(float *)calloc(np,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd,&n);

   y=(float *)calloc(ny,sizeof(float));

   for(int i=0;i<nx;i++)
   {
      for(int l=0;l<ny;l++) y[l]=image[nx*l +i];

      for(int j=0;j<ny;j++)
         image2[nx*j+i]=conv_pnt_sk(y,ny,h,n,j);
   }

   free(h);
   free(y);
   
   return(image2);
}

float *smoothY(float *image, int nx, int ny, int nz, float sd)
{
   int n;
   float *h;
   float *image2;
   float *y;
   int slice_offset;

   image2=(float *)calloc(nx*ny*nz,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image[i];
      return(image2);
   }

   h = gaussian_kernel(sd, &n);
   if(h==NULL)
   {
      return(NULL);
   }

   y=(float *)calloc(ny,sizeof(float));

   for(int k=0;k<nz;k++)
   {
      slice_offset = nx*ny*k;
      for(int i=0;i<nx;i++)
      {
         for(int j=0;j<ny;j++) 
         {
            y[j]=image[slice_offset + nx*j + i];
         }

         for(int j=0;j<ny;j++)
         {
            image2[slice_offset + nx*j + i] = conv_pnt_sk(y,ny,h,n,j);
         }
      }
   }

   free(h);
   free(y);

   return(image2);
}

/////////////////////////////////////////////////////////////////////

float *smoothZ(float *image1, int nx, int ny, int nz, float sd)
{
   int n;
   int np;
   int row_offset;
   float *h;
   float *z;
   float *image2;

   np = nx*ny;

   image2=(float *)calloc(np*nz ,sizeof(float));
   if(image2==NULL) return(NULL);

   if(sd<=0.0)
   {
      for(int i=0; i<nx*ny; i++) image2[i] = image1[i];
      return(image2);
   }

   h = gaussian_kernel(sd,&n);
   if(h==NULL)
   {
      return(NULL);
   }

   z=(float *)calloc(nz,sizeof(float));

   for(int j=0;j<ny;j++)
   {
      row_offset = nx*j;

      for(int i=0;i<nx;i++)
      {
         for(int k=0;k<nz;k++)
         {
            z[k] = image1[np*k + row_offset +i];
         }

         for(int k=0;k<nz;k++)
         {
            image2[np*k + row_offset +i] = conv_pnt_sk(z,nz,h,n,k);
         }
      }
   }

   free(h);
   free(z);

   return(image2);
}

/////////////////////////////////////////////////////////////////////

float *smoothXYZ(float *image1, int nx, int ny, int nz, float sdx, float sdy, float sdz)
{
   float *imA, *imB, *imC;

   imA=smoothX(image1, nx, ny, nz, sdx);
   if(imA==NULL)
   {
      return(NULL);
   }

   imB=smoothY(imA, nx, ny, nz, sdy);
   free(imA);
   if(imB==NULL)
   {
      return(NULL);
   }

   imC=smoothZ(imB, nx, ny, nz, sdz);
   free(imB);

   return(imC);
}

float *smoothXYZ(short *image1, int nx, int ny, int nz, float sdx, float sdy, float sdz)
{
   float *imA, *imB, *imC;

   imA=smoothX(image1, nx, ny, nz, sdx);
   if(imA==NULL)
   {
      return(NULL);
   }

   imB=smoothY(imA, nx, ny, nz, sdy);
   free(imA);
   if(imB==NULL)
   {
      return(NULL);
   }

   imC=smoothZ(imB, nx, ny, nz, sdz);
   free(imB);

   return(imC);
}

/////////////////////////////////////////////////////////////////////
float *smoothXY(float *image1, int nx, int ny, float sd)
{
   float *imA, *imB;

   imA=smoothX(image1, nx, ny, sd);

   imB=smoothY(imA, nx, ny, sd);
   free(imA);

   return(imB);
}

float *smoothXY(short *image1, int nx, int ny, float sd)
{
   float *imA, *imB;

   imA=smoothX(image1, nx, ny, sd);

   imB=smoothY(imA, nx, ny, sd);
   free(imA);

   return(imB);
}
/////////////////////////////////////////////////////////////////////
