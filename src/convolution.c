/**************************************************
   program: convolution.c
   Copyright (c) 1997 Babak A. Ardekani
   ALL RIGHTS RESERVED

   created: 26 June 1997
***************************************************/
#define _convolution

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float conv_pnt_sk(short *x,int sx,float *h,int sh,int i0);
float conv_pnt_sk(float *x,int sx,float *h,int sh,int i0);
float *conv_sk(short *x,int sx,float *h,int sh);
float *conv_sk(float *x,int sx,float *h,int sh);
void conv_sk(float *x, float *y, int sx,float *h,int sh);
void conv_sk_inplace(float *x, int sx, float *h, int sh);

/**************************************************
Function
   float conv_pnt_sk(short *x,int sx,float *h,int sh,int i0)
Inputs 
   x  - an array of type "short" (x[0],x[1],...,x[sx-1])
   sx - size of x 
   h  - an array of type "float" (h[0],h[1],...,h[sh-1])
        containing the right-half of a symmetric (even)
        kernel k (including k[0]), that is:

              k[i] = h[i] i>=0
              k[i] = h[-i] i<0

   sh - size of h (size of k is 2*sh-1)
   i0 - the point at which the convolution sum is evaluated

Outputs
   returned value - convolution of k and x evaluated at i0

Notes
   This function computes and returns \sum_i x[i]*k[i-i0] 
   i.e., the convolution of k and x evaluated at i0.  The
   kernel function k is assumed to be symmetric (even).

   The kernel size 2*sh-1 must be less than or equal to sx, 
   otherwise the function returns 0.

   If i0 is too large/small such that k[i-i0] and x[i] do not overlap
   at some points, array x is extended by reflection in order to 
   compute the convolution sum.
***************************************************/

float conv_pnt_sk(short *x,int sx,float *h,int sh,int i0)
{
	float sum;   /* the convolution sum retured by this function */
	int k;

	if( i0<0 || i0>=sx || (2*sh-1)>sx ) return(0.0);

	sum=0.0;
	if(i0<(sh-1))  
	{
		/* partial overlap: k[i-i0] runs off the left edge of x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
		{
			k=i0-i;
		        if(k<0) k *= (-1);
			sum += h[i]*x[k];
		}
	}
	else if(i0>(sx-sh))
	{
		/* partial overlap: k[i-i0] runs off the right edge of x[i] */
		for(int i=1;i<sh;i++)
		{
			k=i0+i;
			if(k>=sx) k = sx-(k-sx+1);
			sum += h[i]*x[k];
		}

		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0-i];
	}
	else
	{
		/* complete overlap of k[i-i0] and x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
			sum += h[i]*x[i0-i];
	}
	return(sum);
}

float conv_pnt_sk(unsigned char *x,int sx,float *h,int sh,int i0)
{
	float sum;   /* the convolution sum retured by this function */
	int k;

	if( i0<0 || i0>=sx || (2*sh-1)>sx ) return(0.0);

	sum=0.0;
	if(i0<(sh-1))  
	{
		/* partial overlap: k[i-i0] runs off the left edge of x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
		{
			k=i0-i;
		        if(k<0) k *= (-1);
			sum += h[i]*x[k];
		}
	}
	else if(i0>(sx-sh))
	{
		/* partial overlap: k[i-i0] runs off the right edge of x[i] */
		for(int i=1;i<sh;i++)
		{
			k=i0+i;
			if(k>=sx) k = sx-(k-sx+1);
			sum += h[i]*x[k];
		}

		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0-i];
	}
	else
	{
		/* complete overlap of k[i-i0] and x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
			sum += h[i]*x[i0-i];
	}
	return(sum);
}

float conv_pnt_sk(float *x,int sx,float *h,int sh,int i0)
{
	float sum;   /* the convolution sum retured by this function */
	int k;

	if( i0<0 || i0>=sx || (2*sh-1)>sx ) return(0.0);

	sum=0.0;
	if(i0<(sh-1))  
	{
		/* partial overlap: k[i-i0] runs off the left edge of x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
		{
			k=i0-i;
			if(k<0) k *= (-1);
			sum += h[i]*x[k];
		}
	}
	else if(i0>(sx-sh))
	{
		/* partial overlap: k[i-i0] runs off the right edge of x[i] */
		for(int i=1;i<sh;i++)
		{
			k=i0+i;
			if(k>=sx) k = sx-(k-sx+1);
			sum += h[i]*x[k];
		}

		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0-i];
	}
	else
	{
		/* complete overlap of k[i-i0] and x[i] */
		for(int i=0;i<sh;i++)
			sum += h[i]*x[i0+i];

		for(int i=1;i<sh;i++)
			sum += h[i]*x[i0-i];
	}

	return(sum);
}

/**************************************************
Function
   float *conv_sk(short *x,int sx,float *h,int sh)
Inputs
   x  - an array of type "short" (x[0],x[1],...,x[sx-1])
   sx - size of x
   h  - an array of type "float" (h[0],h[1],...,h[sh-1])
        containing the right-half of a symmetric (even)
        kernel k (including k[0]), that is:

              k[i] = h[i] i>=0
              k[i] = h[-i] i<0

   sh - size of h (size of k is 2*sh-1)

Outputs
   returned array - convolution of k and x 

Notes
   This function computes and returns \sum_i x[i]*k[i-i0]
   for i0=0,1,...,sx, (i.e., convolution of x and k). 
   The kernel function k is assumed to be symmetric (even).

   The kernel size 2*sh-1 must be less than or equal to sx,
   otherwise the function returns 0's.

   Returns NULL in case of memory allocation problem.
***************************************************/

float *conv_sk(short *x,int sx,float *h,int sh)
{
	int i;
	float *y;

	y=(float *)calloc(sx,sizeof(float)); 
	if(y==NULL) return(NULL);

	for(i=0;i<sx;i++)
		y[i]=conv_pnt_sk(x,sx,h,sh,i);

	return(y);
}

float *conv_sk(float *x,int sx,float *h,int sh)
{
	int i;
	float *y;

	y=(float *)calloc(sx,sizeof(float)); 
	if(y==NULL) return(NULL);

	for(i=0;i<sx;i++)
		y[i]=conv_pnt_sk(x,sx,h,sh,i);

	return(y);
}

void conv_sk(float *x, float *y, int sx,float *h,int sh)
{
	for(int i=0;i<sx;i++)
		y[i]=conv_pnt_sk(x,sx,h,sh,i);
}

void conv_sk_inplace(float *x, int sx, float *h, int sh)
{
	int i;
	float *y;

	y=(float *)calloc(sx,sizeof(float)); 
	if(y==NULL) return;

	for(i=0;i<sx;i++)
		y[i]=conv_pnt_sk(x,sx,h,sh,i);

	for(i=0;i<sx;i++) x[i]=y[i];

	free(y);
}
