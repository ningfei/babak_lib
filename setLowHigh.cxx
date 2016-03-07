#include <stdlib.h>
#include "minmax.h"

void setLowHigh(short *image, int nv, int *low, int *high)
{
   short min, max;
   int *histogram;
   int hsize;			/* histogram size */
   int b;
   int i;

   int nmax;
   int n;

   minmax(image,nv,min,max);

   hsize = max-min+1;
   
   histogram=(int *)calloc(hsize,sizeof(int));

   for(i=0;i<nv;i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int)(0.0005 * nv);

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[i];
      if(n>nmax) break;
   }
   *low=i+min; 

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[hsize-1-i];
      if(n>nmax) break;
   }
   *high=hsize-1-i+min; 

   free(histogram);
}

void setLowHigh(short *image, int nv, int *low, int *high, float percent)
{
   short min, max;
   int *histogram;
   int hsize;			/* histogram size */
   int b;
   int i;

   int nmax;
   int n;

   minmax(image,nv,min,max);

   hsize = max-min+1;

   histogram=(int *)calloc(hsize,sizeof(int));

   for(int i=0; i<nv; i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int)( percent * nv/100.0);

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[i];
      if(n>nmax) break;
   }
   *low=i+min; 

   n=0;
   for(i=0;i<hsize;i++)
   {
      n += histogram[hsize-1-i];
      if(n>nmax) break;
   }
   *high=hsize-1-i+min; 

   free(histogram);
}
