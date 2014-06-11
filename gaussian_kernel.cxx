/**************************************************
   program: gaussian_kernel.c
   Copyright (c) 1997 Babak A. Ardekani
   ALL RIGHTS RESERVED

   created: 25 June 1997
   revision: 26 June 1997 (BAA)
             gaussian_kernel() modified to handle the case were the intput sd<=0 
   revision: 3 August 2010 (BAA)
             removed the third input parameter (v) for printing
***************************************************/

#define _gaussian_kernel

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**************************************************
Inputs 
   sd - standard deviation of the Gaussian IN UNITS OF PIXELS 

Outputs
   n - dimension of the returned array
   returned pointer - array of dimension n (h[0],h[1],...,h[n-1])

Notes
   This function computes and returns an array of values proportional to:
   exp( -0.5*x^2/sd^2 ) ( x=0, 1, 2, ..., n-1 ).

   Since Gaussian curves are symmetric, only the right-half of the curve 
   (including 0) is computed and returned (n points).  The complete curve 
   would have 2*n-1 points.  

   The returned array is scaled such that the sum of all the coefficients 
   (including the left-half of the curve) would be 1.

   The dimension n is chosen such that the area under the Gaussian curve is 
   at least 0.99 
***************************************************/

float *gaussian_kernel(float sd, int *n)
{
   float *h;  /* the array to be returned */
   float var; /* constant variance */
   float nf;  /* normalization factor */

   // sd cannot be negative
   if(sd<0.0)
   {
      sd = 0.0;
      *n=1;
      h=(float *)calloc(*n,sizeof(float));
      h[0]=1.0;
      return(h);
   }

   var=sd*sd;

   // factor 2.57 was chosen to make the area under the Gaussian at least 0.99 
   // factor 1.65 would give an area of 0.90
   *n=(int)ceilf(sd*2.57);

   // ensure that n is not zero
   if(*n==0) 
   {
      *n=1;	
   }

   h=(float *)calloc(*n,sizeof(float));
   if(h==NULL) 
   {
      return(NULL);
   }
  
   nf=0.0;
   for(int x=0; x<(*n); x++)
   {
      h[x] = expf(-0.5*x*x/var);
      nf += h[x];
   }

   /* add the values of the left-half of the Gaussian curve to nf */
   nf = 2.0*nf - h[0];

   /* normalize the Gaussian curve to 1 */
   for(int x=0; x<(*n); x++)
   {
      h[x] /= nf;
   }

   return(h);
}
