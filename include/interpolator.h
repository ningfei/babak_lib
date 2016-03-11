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


template<class TYPE> float linearInterpolator(float x, float y, float z, TYPE *array, int nx, int ny, int nz, int np)
{
   int n_and_nx;  // n + nx
   int n_and_np;  // n + np;

   int   i,j,k,n;
   float u1,u2,u3,u4,u5,u6,u7,u8;
   float v1,v2,v3,v4;
   float w1,w2;
   float xr;
   float yr;
   float zr;

   i=(int)(x);
   j=(int)(y);
   k=(int)(z);

   if(i<0 || i>(nx-2) || j<0 || j>(ny-2) )
   {
      return(0.0);
   }

   if( k>=0 && k<(nz-1) )
   {
      xr = x - i;

      n= k*np + j*nx + i;

      n_and_nx = n + nx;
      n_and_np = n + np;

      u1 = array[n]; 
      u2 = array[n+1]; 
      u3 = array[n_and_nx];
      u4 = array[n_and_nx + 1];
      u5 = array[n_and_np];
      u6 = array[n_and_np + 1];
      u7 = array[n_and_nx + np];
      u8 = array[n_and_nx + np + 1];

      v1 = u1 + (u2-u1)*xr;
      v2 = u3 + (u4-u3)*xr;
      v3 = u5 + (u6-u5)*xr;
      v4 = u7 + (u8-u7)*xr;

      yr = y - j;
      w1 = v1 + (v2-v1)*yr;
      w2 = v3 + (v4-v3)*yr;

      return( w1 + (w2-w1)*(z-k) ); // saved a multiplication :)
   }
   else if( k==(nz-1) )
   {
      xr = x - i;

      n=k*np + j*nx +i;

      v1 = array[n];
      v2 = array[n+1];
      w1 = v1 + (v2-v1)*xr;

      n += nx;

      v1 = array[n];
      v2 = array[n+1];

      w2 = v1 + (v2-v1)*xr;

      return( w1 + (w2-w1)*(y-j) );
   }

   return(0.0);
}

template<class TYPE> float nearestNeighbor(float x, float y, float z, TYPE *array, int nx, int ny, int nz, int np)
{
   int i,j,k;
  
   i=(int)(x+0.5);
   j=(int)(y+0.5);
   k=(int)(z+0.5);

   if( i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz)
   {
      return( (float)(array[ np*k + nx*j +i ]) );
   }
   else
   {
      return(0.0);
   }
}

#define _linearInterpolator_h
#endif
