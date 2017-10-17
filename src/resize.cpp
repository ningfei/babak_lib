#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <babak_lib.h>
#include <resizeX.h>

int2 *resizeZ(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dz1, int4 nz2, float4 dz2);
float4 *resizeZ(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dz1, int4 *nz2, float4 *dz2);

float4 *resizeY(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dy1, int4 ny2, float4 dy2);
int2 *resizeY(float4 *image1, int4 nx1, int4 ny1, float4 dy1, int4 ny2, float4 dy2);
float4 *resizeY(float4 *image1, int4 nx1, int4 ny1, float4 dy1, int4 ny2, float4 *dy2);

float4 *resizeXY(float4 *image1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2)
{
	float4 *imA;
	float4 *image2;

	imA = resizeX(image1, nx1, ny1, dx1, nx2, dx2);

	image2 = resizeY(imA, nx2, ny1, dy1, ny2, &dy2);
	free(imA);
   
	return(image2);
}

float4 *resizeY(float4 *image1, int4 nx1, int4 ny1, float4 dy1, int4 ny2, float4 *dy2)
{
   int4 n;
   int4 np2;
   int4 j0;
   float4 d,dl;
   float4 sd=0.0;
   float4 *h;
   float4 *image2;
   float4 *x;
   float4 yc1,yc2;

   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = (*dy2)*(ny2-1.0)/2.0;

   image2=(float4 *)calloc(np2,sizeof(float4));
   if(image2==NULL) return(NULL);

   if(dy1 < (*dy2) )
   {
      sd=(float4)sqrt( (0.5/log(2.0))*( (*dy2)*(*dy2) - dy1*dy1 ) );
      sd /= dy1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float4 *)calloc(ny1,sizeof(float4));

   for(int4 i=0;i<nx1;i++)
   {
      for(int4 l=0;l<ny1;l++)
         x[l]=image1[nx1*l +i];

      for(int4 j=0;j<ny2;j++)
      {
         d = (j*(*dy2) - yc2 + yc1)/dy1;
         j0=(int4)d;
         dl=d - j0;

         if(dl==0.0)
            image2[nx1*j+i]=(conv_pnt_sk(x,ny1,h,n,j0));
         else
            image2[nx1*j+i]=( (1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1) );
      }
   }

   free(h);
   free(x);

   return(image2);
}

int2 *resizeXY(int2 *image1, int4 nx1, int4 ny1, float4 dx1, float4 dy1, int4 nx2, int4 ny2, float4 dx2, float4 dy2)
{
	float4 *imA;
	int2 *image2;

	imA = resizeX(image1, nx1, ny1, dx1, nx2, dx2);

	image2 = resizeY(imA, nx2, ny1, dy1, ny2, dy2);
	free(imA);
   
	return(image2);
}

int2 *resizeY(float4 *image1, int4 nx1, int4 ny1, float4 dy1, int4 ny2, float4 dy2)
{
   int4 n;
   int4 np2;
   int4 j0;
   float4 d,dl;
   float4 sd=0.0;
   float4 *h;
   int2 *image2;
   float4 *x;
   float4 yc1,yc2;

   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = dy2*(ny2-1.0)/2.0;

   image2=(int2 *)calloc(np2,sizeof(int2));
   if(image2==NULL) return(NULL);

   if(dy1 < dy2)
   {
      sd=(float4)sqrt( (0.5/log(2.0))*( dy2*dy2 - dy1*dy1 ) );
      sd /= dy1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float4 *)calloc(ny1,sizeof(float4));

   for(int4 i=0;i<nx1;i++)
   {
      for(int4 l=0;l<ny1;l++)
         x[l]=image1[nx1*l +i];

      for(int4 j=0;j<ny2;j++)
      {
         d = (j*dy2 - yc2 + yc1)/dy1;
         j0=(int4)d;
         dl=d - j0;

         if(dl==0.0)
            image2[nx1*j+i]=(int2)(conv_pnt_sk(x,ny1,h,n,j0) + 0.5);
         else
            image2[nx1*j+i]=(int2)( (1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1) + 0.5 );
      }
   }

   free(h);
   free(x);

   return(image2);
}

int2 *resizeXYZ(char *image1, 
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2)
{
	float4 *imA, *imB;
	int2 *image2;
	int2 *tmp;
	int4 nv1;

	nv1 = nx1*ny1*nz1;

	tmp = (int2 *)calloc(nv1,sizeof(int2));
	for(int4 i=0; i<nv1; i++) tmp[i] = image1[i];
	imA = resizeX(tmp, nx1, ny1, nz1, dx1, nx2, dx2);
	free(tmp);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

uchar *resizeZ_UC(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dz1, int4 nz2, float4 dz2)
{
	int4 n;
	int4 np1;
	int4 k0;
	float4 d,dl;
	float4 sd=0.0;
	float4 *h;
	float4 *x;
	float4 zc1,zc2;
	uchar *image2;

	np1=nx1*ny1;

	zc1 = dz1*(nz1-1.0)/2.0;
	zc2 = dz2*(nz2-1.0)/2.0;

	image2=(uchar *)calloc(np1*nz2,1);
	if(image2==NULL) return(NULL);

	if(nz2 < nz1)
	{
		sd=(float4)sqrt( (0.5/log(2.0))*( dz2*dz2 - dz1*dz1 ) );

		sd /= dz1;

		h = gaussian_kernel(sd,&n);

		x=(float4 *)calloc(nz1,sizeof(float4));

		for(int4 j=0;j<ny1;j++)
		for(int4 i=0;i<nx1;i++)
		{
			for(int4 l=0;l<nz1;l++)
				x[l]=image1[np1*l + nx1*j +i];

			for(int4 k=0;k<nz2;k++)
			{
				d = (k*dz2 - zc2 + zc1)/dz1;
				k0=(int4)d;
				dl=d - k0;

				if(dl==0)
					image2[np1*k+nx1*j+i]=(uchar)(conv_pnt_sk(x,nz1,h,n,k0)+0.5);
				else
					image2[np1*k+nx1*j+i]=(uchar) ( (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) +
                   	dl*conv_pnt_sk(x,nz1,h,n,k0+1) + 0.5 );
			}
		}

		free(h);
		free(x);
	}
	else if(nz2 > nz1)
	{
		for(int4 j=0;j<ny1;j++)
		for(int4 i=0;i<nx1;i++)
		for(int4 k=0;k<nz2;k++)
		{
			d = (k*dz2 - zc2 + zc1)/dz1;
			k0=(int4)d;
			dl=d - k0;

			if(k0<0 || k0>=(nz1-1))
			{
				image2[np1*k+nx1*j+i]=0;
				continue;
			}

			if(dl==0.0)
				image2[np1*k+nx1*j+i]=(uchar)(image1[np1*k0 + nx1*j +i]+0.5);
			else
				image2[np1*k+nx1*j+i]=(uchar)((1.0-dl)*image1[np1*k0 + nx1*j +i]
					+ dl*image1[np1*(k0+1) + nx1*j +i] + 0.5 );
		}
	}
	else
	{
		for(int4 i=0;i<np1*nz1;i++)
			image2[i]=(uchar)(image1[i]+0.5);
	}

   return(image2);
}

uchar *resizeXYZ(uchar *image1, 
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2)
{
	float4 *imA, *imB;
	uchar *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ_UC(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

int2 *resizeXYZ(int2 *image1, 
int4 nx1, int4 ny1, int4 nz1, float4 dx1, float4 dy1, float4 dz1,
int4 nx2, int4 ny2, int4 nz2, float4 dx2, float4 dy2, float4 dz2)
{
	float4 *imA, *imB;
	int2 *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, nz2, dz2);
	free(imB);

	return(image2);
}

int2 *resizeXYZ(int2 *image1,  DIM dim1, DIM dim2)
{
	float4 *imA, *imB;
	int2 *image2;

	imA = resizeX(image1, dim1.nx, dim1.ny, dim1.nz, dim1.dx, dim2.nx, dim2.dx);

	imB = resizeY(imA, dim2.nx, dim1.ny, dim1.nz, dim1.dy, dim2.ny, dim2.dy);
	free(imA);
   
	image2=resizeZ(imB, dim2.nx, dim2.ny, dim1.nz, dim1.dz, dim2.nz, dim2.dz);
	free(imB);

	return(image2);
}

int2 *resizeZ(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dz1, int4 nz2, float4 dz2)
{
   int4 n;
   int4 np1;
   int4 k0;
   float4 d,dl;
   float4 sd=0.0;
   float4 *h;
   float4 *x;
   float4 zc1,zc2;
   int2 *image2;

   np1=nx1*ny1;

   zc1 = dz1*(nz1-1.0)/2.0;
   zc2 = dz2*(nz2-1.0)/2.0;

   image2=(int2 *)calloc(np1*nz2,sizeof(int2));
   if(image2==NULL) return(NULL);

   if(dz1 < dz2)
   {
      sd=(float4)sqrt( (0.5/log(2.0))*( dz2*dz2 - dz1*dz1 ) );
      sd /= dz1;
   }
   else sd = 0.0;

   h = gaussian_kernel(sd,&n);

   x=(float4 *)calloc(nz1,sizeof(float4));

   for(int4 j=0;j<ny1;j++)
   for(int4 i=0;i<nx1;i++)
   {
      for(int4 l=0;l<nz1;l++)
         x[l]=image1[np1*l + nx1*j +i];

      for(int4 k=0;k<nz2;k++)
      {
         d = (k*dz2 - zc2 + zc1)/dz1;
         k0=(int4)d;
         dl=d - k0;
         
         if(dl==0)
            image2[np1*k+nx1*j+i]=(int2)(conv_pnt_sk(x,nz1,h,n,k0)+0.5);
         else
            image2[np1*k+nx1*j+i]=(int2) ( (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) +
            dl*conv_pnt_sk(x,nz1,h,n,k0+1) + 0.5 );
      }
   }

   free(h);
   free(x);

   return(image2);
}

float4 *resizeXYZ(float4 *image1, 
	int4 nx1, int4 ny1, int4 nz1, 
	float4 dx1, float4 dy1, float4 dz1,
	int4 nx2, int4 ny2, int4 nz2, 
	float4 dx2, float4 dy2, float4 dz2)
{
	float4 *imA, *imB;
	float4 *image2;

	imA = resizeX(image1, nx1, ny1, nz1, dx1, nx2, dx2);

	imB = resizeY(imA, nx2, ny1, nz1, dy1, ny2, dy2);
	free(imA);
   
	image2=resizeZ(imB, nx2, ny2, nz1, dz1, &nz2, &dz2);
	free(imB);

	return(image2);
}

float4 *resizeY(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dy1, int4 ny2, float4 dy2)
{
   int4 n;
   int4 np1, np2;
   int4 j0;
   float4 d,dl;
   float4 sd=0.0;
   float4 *h;
   float4 *image2;
   float4 *x;
   float4 yc1,yc2;

   np1=nx1*ny1;
   np2=nx1*ny2;

   yc1 = dy1*(ny1-1.0)/2.0;
   yc2 = dy2*(ny2-1.0)/2.0;

   image2=(float4 *)calloc(np2*nz1,sizeof(float4));
   if(image2==NULL) return(NULL);

   if(dy1 < dy2)
   {
      sd=(float4)sqrt( (0.5/log(2.0))*( dy2*dy2 - dy1*dy1 ) );
      sd /= dy1;
   }
   else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float4 *)calloc(ny1,sizeof(float4));

   for(int4 k=0;k<nz1;k++)
   for(int4 i=0;i<nx1;i++)
   {
      for(int4 l=0;l<ny1;l++)
         x[l]=image1[np1*k + nx1*l +i];

      for(int4 j=0;j<ny2;j++)
      {
         d = (j*dy2 - yc2 + yc1)/dy1;
         j0=(int4)d;
         dl=d - j0;

         if(dl==0.0)
            image2[np2*k+nx1*j+i]=conv_pnt_sk(x,ny1,h,n,j0);
         else
            image2[np2*k+nx1*j+i]=(1.0-dl)*conv_pnt_sk(x,ny1,h,n,j0) + dl*conv_pnt_sk(x,ny1,h,n,j0+1);
      }
   }

   free(h);
   free(x);

   return(image2);
}

float4 *resizeZ(float4 *image1, int4 nx1, int4 ny1, int4 nz1, float4 dz1, int4 *nz2, float4 *dz2)
{
   int4 n;
   int4 np1;
   int4 k0;
   float4 d,dl;
   float4 sd=0.0;
   float4 *h;
   float4 *x;
   float4 zc1,zc2;
   float4 *image2;

   np1=nx1*ny1;

   zc1 = dz1*(nz1-1.0)/2.0;
   zc2 = (*dz2)*( (*nz2) -1.0)/2.0;

   image2=(float4 *)calloc(np1*(*nz2) ,sizeof(float4));
   if(image2==NULL) return(NULL);

   if( dz1 < (*dz2) )
   {
      sd=(float4)sqrt( (0.5/log(2.0))*( (*dz2)*(*dz2) - dz1*dz1 ) );
      sd /= dz1;
   } else sd=0.0;

   h = gaussian_kernel(sd,&n);

   x=(float4 *)calloc(nz1,sizeof(float4));

   for(int4 j=0;j<ny1;j++)
   for(int4 i=0;i<nx1;i++)
   {
      for(int4 l=0;l<nz1;l++)
         x[l]=image1[np1*l + nx1*j +i];

      for(int4 k=0;k<(*nz2);k++)
      {
         d = (k*(*dz2) - zc2 + zc1)/dz1;
         k0=(int4)d;
         dl=d - k0;

         if(dl==0)
            image2[np1*k+nx1*j+i]=conv_pnt_sk(x,nz1,h,n,k0);
         else
            image2[np1*k+nx1*j+i]= (1.0-dl)*conv_pnt_sk(x,nz1,h,n,k0) + dl*conv_pnt_sk(x,nz1,h,n,k0+1);
      }
   }

   free(h);
   free(x);

   return(image2);
}
