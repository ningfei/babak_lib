#include <babak_lib.h>

void compute_cm(short *image, int nx, int ny, int nz, float dx, float dy, float dz, float *x, float *y, float *z)
{
   double dum;
   int q;

   *x=0.0; *y=0.0; *z=0.0;

   dum=0.0;
   q=0;
   for(int k=0;k<nz;k++)
   for(int j=0;j<ny;j++)
   for(int i=0;i<nx;i++)
   {
      dum += image[q];
      *x += i*image[q];
      *y += j*image[q];
      *z += k*image[q];
      q++;
   }

   *x = *x/dum;
   *y = *y/dum;
   *z = *z/dum;

/* transfer origin to the volume center */
   *x -= (nx-1.0)/2.0;
   *y -= (ny-1.0)/2.0;
   *z -= (nz-1.0)/2.0;

/* convert to mm */
   *x *= dx;
   *y *= dy;
   *z *= dz;
}
