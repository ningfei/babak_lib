#ifndef _linearInterpolator_h
#include <cstddef>

template<class TYPE> float linearInterpolator(float x, float y, TYPE *array, int nx, int ny)
{
   int i,j,n;
   float u,uu;
   float v1, v2;
	
   i=(int)(x);
   j=(int)(y);

   if(i<0 || i>(nx-2) || j<0 || j>(ny-2) ) 
   {
      return(0.0);
   }

   u = x - i; if(u<0.0) u=0.0;
   uu = 1.0-u;

   n= j*nx +i;
   v1 = array[n]*uu + array[n+1]*u;
   v2 = array[n+nx]*uu + array[n+nx+1]*u;

   u = y - j; if(u<0.0) u=0.0;
   uu = 1.0-u;

   return( v1*uu + v2*u );
}
#define _linearInterpolator_h
#endif
