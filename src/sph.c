#include "../include/babak_lib.h"
#include "../include/sph.h"
#include <math.h>
#include <stdlib.h>

float4 SPH::dot(float4 *u)
{
   float4 dot = 0.0;

   for(int4 m=0; m<n; m++)
      dot += (v[m]*u[m]);

   return(dot);
}

float4 SPH::norm()
{
   float4 norm=0.0;

   for(int4 m=0; m<n; m++)
      norm += (v[m]*v[m]);

   norm = sqrtf(norm);

   return(norm);
}

float4 SPH::mean()
{
   float4 mean=0.0;

   for(int4 m=0; m<n; m++)
      mean += v[m];

   mean /= n;

   return(mean);
}

void SPH::normalize()
{
   float4 norm=0.0;

   for(int4 m=0; m<n; m++)
      norm += (v[m]*v[m]);

   norm = sqrtf(norm);

   for(int4 m=0; m<n; m++)
      v[m] /= norm;
}

void SPH::zeromean()
{
   float4 mean=0.0;

   for(int4 m=0; m<n; m++)
      mean += v[m];

   mean /= n;

   for(int4 m=0; m<n; m++)
      v[m] -= mean;
}

void SPH::get(float4 *array)
{
   for(int4 m=0; m<n; m++)
   {
      array[m] = v[m];
   }
}

void SPH::set(SHORTIM im, int4 ic, int4 jc, int4 kc)
{
   int4 i0, j0, k0;

   for(int4 m=0; m<n; m++)
   {
      i0=ic+i[m];
      j0=jc+j[m];
      k0=kc+k[m];

      if(i0>=0 && i0<im.nx && j0>=0 && j0<im.ny && k0>=0 && k0<im.nz)
      {
         v[m] = im.v[ im.np*k0 + im.nx*j0 + i0 ];
      } else
      {
         v[m] = 0.0;
      }
   }
}

void SPH::reset()
{
   for(int4 m=0; m<n; m++)
   {
      v[m]=0.0;
   }
}

SPH::~SPH()
{
   delete i;
   delete j;  
   delete k;  
   delete v;  
}

SPH::SPH(int4 radius)
{
   r=radius;

   int4 nmax;
   float4 r2;
   float4 x2, y2;

   if(r<0) r=0;

   nmax = (2*r+1) * (2*r+1) * (2*r+1);

   i = (int4 *)calloc(nmax, sizeof(int4));
   j = (int4 *)calloc(nmax, sizeof(int4));
   k = (int4 *)calloc(nmax, sizeof(int4));

   r2 = r*r;

   n=0;
   for(int4 x=-r; x<=r; x++)   
   {
      x2=x*x;
      for(int4 y=-r; y<=r; y++)   
      {
         y2=y*y;
         for(int4 z=-r; z<=r; z++)   
         {
            if( (x2+y2+z*z) <= r2 )
            {
               i[n]=x; j[n]=y; k[n]=z;
               n++;
            }
         }
      }
   }

   v = (float4 *)calloc(n, sizeof(float4));
}

