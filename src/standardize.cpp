#include <math.h>

void standardize(float *x, int n)
{
   float v;
   double mean;
   double norm;
   double s=0.0;
   double s2=0.0;

   for(int i=0; i<n; i++)
   {
      v = x[i];
      s += v;
      s2 += v*v;
   }

   mean = s/n;
   norm = sqrt( s2 - s*s/n );

   for(int i=0; i<n; i++)
   {
      x[i] = (float) ( (x[i] - mean)/norm );
   }
}

void standardize(float *x, float *mask, int n)
{
   float v;
   double mean=0.0;
   double norm=1.0;
   double s=0.0;
   double s2=0.0;
   int N=0;

   for(int i=0; i<n; i++)
   if(mask[i]!=0.0)
   {
      v = x[i];
      s += v;
      s2 += v*v;
      N++;
   }

   if(N>0)
   {
      mean = s/N;
      norm = sqrt( s2 - s*s/N );
   }
   else
   {
      mean=0.0;
      norm=1.0;
   }

   for(int i=0; i<n; i++)
   if(mask[i]!=0.0)
   {
      x[i] = (float) ( (x[i] - mean)/norm );
   }
   else
   {
      x[i] = 0.0;
   }
}
