// Revision:
// 2/4/2011 - added slice_seek and row_seek variables to reduce the number of operations
// in the tripple loop.
// Removed the custom hpsort function in this file.
// Modified medianFilter so that the hpsort function in hpsort.c is used, 
// which indexes from 0 rather than 1.
// Added the nv variable.
// Added memory allocation check.

#define _medianfilter

#include <stdlib.h>
#include <stdio.h>
#include "../include/babak_lib.h"

// Note: Wx, Wy, and Wz can be 0
// Replaces image_in with it's median filtered version
void medianFilter(float *image_in, int nx, int ny, int nz, int Wx, int Wy, int Wz)
{
   int np, nv;
   np = nx*ny;
   nv = np*nz;

   // N is always an odd number because (2*n+1) is odd for all n, and the product of
   // odd numbers is always odd.
   int N; // dimension of array (size of 3D window)
   N = (2*Wx+1)*(2*Wy+1)*(2*Wz+1);

   float *array;
   array = (float *)calloc(N,sizeof(float));

   float *image_out;
   image_out = (float *)calloc(nv,sizeof(float));

   if(array==NULL || image_out==NULL)
   {
      printf("\nWarning: memory allocation error in medianFilter()\n");
      return;
   }

   // silly me, I had M equal to the following monster, it is simply (N-1)/2
   // M = (2*Wx+1)*(2*Wy+1)*Wz + 2*Wx*Wy + Wx + Wy;
   // Since N is odd, (N-1) is even (divisible by 2 with no remainders)
   int M; // points to the middle index of the array (i.e., (N-1)/2 )
   M = (N-1)/2;

   int slice_seek, row_seek;
   for(int k=0; k<nz; k++)
   {
      slice_seek = k*np;

      for(int j=0; j<ny; j++)
      {
         row_seek = slice_seek + j*nx;

         for(int i=0; i<nx; i++)
         {
            extractArray(image_in, nx, ny, nz, np, i, j, k, Wx, Wy, Wz, array);
            hpsort(N,array);
            image_out[row_seek + i]=array[M];
         }
      }
   }
   
   for(int v=0; v<nv; v++)
   {
      image_in[v]=image_out[v];
   }

   free(image_out);
}
