#include <stdlib.h>
#include "../include/minmax.h"
#include "../include/babak_lib.h"

void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high)
{
   int2 min, max;
   int4 *histogram;
   int4 hsize;			/* histogram size */
   int4 b;
   int4 i;

   int4 nmax;
   int4 n;

   minmax(image,nv,min,max);

   hsize = max-min+1;
   
   histogram=(int4 *)calloc(hsize,sizeof(int4));

   for(i=0;i<nv;i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int4)(0.0005 * nv);

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

void setLowHigh(int2 *image, int4 nv, int4 *low, int4 *high, float4 percent)
{
   int2 min, max;
   int4 *histogram;
   int4 hsize;			/* histogram size */
   int4 b;
   int4 i;

   int4 nmax;
   int4 n;

   minmax(image,nv,min,max);

   hsize = max-min+1;

   histogram=(int4 *)calloc(hsize,sizeof(int4));

   for(int4 i=0; i<nv; i++)
   {
      b = image[i] - min;

      if(b>=0 && b<hsize) histogram[ b ]++;
   }

   nmax = (int4)( percent * nv/100.0);

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
