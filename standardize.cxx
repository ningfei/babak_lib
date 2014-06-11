#include <babak_lib.h>

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
