#include <babak_lib.h>

short *reorientVolume(short *v1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1, float *orientMat,
int *nx2, int *ny2, int *nz2, float *dx2, float *dy2, float *dz2)
{
   short *v2;
   int np1, nv1;
   int np2, nv2;
   int i1_a, j1_a, k1_a; 
   int i1_b, j1_b, k1_b; 
   int i1, j1, k1; 
   int xc,yc,zc;

   np1 = nx1*ny1;
   nv1 = nx1*ny1*nz1;

   v2  = (short *)calloc(nv1, sizeof(short));
   if(v2==NULL)
   {
      memory_allocation_error("v2");
   }

   if     ( orientMat[0] != 0.0 ) 	{ *nx2 = nx1; *dx2 = dx1; }
   else if( orientMat[4] != 0.0 ) 	{ *ny2 = nx1; *dy2 = dx1; }
   else if( orientMat[8] != 0.0 ) 	{ *nz2 = nx1; *dz2 = dx1; }

   if     ( orientMat[1] != 0.0 ) 	{ *nx2 = ny1; *dx2 = dy1; }
   else if( orientMat[5] != 0.0 ) 	{ *ny2 = ny1; *dy2 = dy1; }
   else if( orientMat[9] != 0.0 ) 	{ *nz2 = ny1; *dz2 = dy1; }

   if     ( orientMat[2] != 0.0 ) 	{ *nx2 = nz1; *dx2 = dz1; }
   else if( orientMat[6] != 0.0 ) 	{ *ny2 = nz1; *dy2 = dz1; }
   else if( orientMat[10] != 0.0)  { *nz2 = nz1; *dz2 = dz1; }

   np2 = (*nx2) * (*ny2);
   nv2 = (*nx2) * (*ny2) * (*nz2);

   if ( orientMat[0] < 0.0 || orientMat[4] < 0.0 || orientMat[8] < 0.0 )
      xc=nx1-1;
   else
      xc=0;

   if ( orientMat[1] < 0.0 || orientMat[5] < 0.0 || orientMat[9] < 0.0 )
      yc=ny1-1;
   else
      yc=0;

   if ( orientMat[2] < 0.0 || orientMat[6] < 0.0 || orientMat[10] < 0.0 )
      zc=nz1-1;
   else
      zc=0;

   int q=0;
   for(int k2=0; k2<(*nz2); k2++)
   {
      i1_a = (int)orientMat[8]*k2 + xc;
      j1_a = (int)orientMat[9]*k2 + yc;
      k1_a = (int)orientMat[10]*k2 + zc;
      for(int j2=0; j2<(*ny2); j2++)
      {
         i1_b = i1_a + (int)orientMat[4]*j2;
         j1_b = j1_a + (int)orientMat[5]*j2;
         k1_b = k1_a + (int)orientMat[6]*j2;
         for(int i2=0; i2<(*nx2); i2++)
         {
            i1 = (int)orientMat[0]*i2 + i1_b;
            j1 = (int)orientMat[1]*i2 + j1_b;
            k1 = (int)orientMat[2]*i2 + k1_b;

            v2[q++]= v1[k1*np1 + j1*nx1 + i1];
         }
      }
   }

   return(v2);
}
