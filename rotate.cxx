#include <babak_lib.h>

/* produces a 4x4 transformation matrix for rotating a point by an angle
alpha about the (x,y,z) axis. */

void rotate(float *R, float alpha, float x, float y, float z)
{
   float CosAlpha,SinAlpha;
   float dum;
   float ax,ay,az;

   R[0]=R[5]=R[10]=R[15]=1.0;
   R[1]=R[2]=R[3]=0.0;
   R[4]=R[6]=R[7]=0.0;
   R[8]=R[9]=R[11]=0.0;
   R[12]=R[13]=R[14]=0.0;
   
   if(x==0.0 && y==0.0 && z==0.0) return;

   CosAlpha=(float)cos((double)alpha);
   SinAlpha=(float)sin((double)alpha);

   /* find the unit vector (ax,ay,az) in the direction of (x,y,z) */
   dum=(float)sqrt((double)x*x+y*y+z*z);
   ax=x/dum; ay=y/dum; az=z/dum;
      
   R[0] = ax*ax + CosAlpha - CosAlpha*ax*ax; 
   R[1] = ax*ay - CosAlpha*ax*ay - SinAlpha*az; 
   R[2] = ax*az - CosAlpha*ax*az + SinAlpha*ay; 
   
   R[4] = ay*ax - CosAlpha*ay*ax + SinAlpha*az; 
   R[5] = ay*ay + CosAlpha - CosAlpha*ay*ay; 
   R[6] = ay*az - CosAlpha*ay*az - SinAlpha*ax; 
   
   R[8]  = az*ax - CosAlpha*az*ax - SinAlpha*ay; 
   R[9]  = az*ay - CosAlpha*az*ay + SinAlpha*ax; 
   R[10] = az*az + CosAlpha - CosAlpha*az*az; 
   
   return;
}

float *rotate(float alpha, float x, float y, float z)
{
   float CosAlpha,SinAlpha;
   float dum;
   float ax,ay,az;
   float *T;

   T=(float *)calloc(16,sizeof(float));
   T[15]=1.0;
   
   if(x==0.0 && y==0.0 && z==0.0)
   {
      T[0]=T[5]=T[10]=1.0;    
      return(T);
   }

   CosAlpha=(float)cos((double)alpha);
   SinAlpha=(float)sin((double)alpha);

   /* find the unit vector (ax,ay,az) in the direction of (x,y,z) */
   dum=(float)sqrt((double)x*x+y*y+z*z);
   ax=x/dum; ay=y/dum; az=z/dum;
   
   T[0] = ax*ax + CosAlpha - CosAlpha*ax*ax; 
   T[1] = ax*ay - CosAlpha*ax*ay - SinAlpha*az; 
   T[2] = ax*az - CosAlpha*ax*az + SinAlpha*ay; 

   T[4] = ay*ax - CosAlpha*ay*ax + SinAlpha*az; 
   T[5] = ay*ay + CosAlpha - CosAlpha*ay*ay; 
   T[6] = ay*az - CosAlpha*ay*az - SinAlpha*ax; 

   T[8]  = az*ax - CosAlpha*az*ax - SinAlpha*ay; 
   T[9]  = az*ay - CosAlpha*az*ay + SinAlpha*ax; 
   T[10] = az*az + CosAlpha - CosAlpha*az*az; 

   return(T);
}

void rotate(float *R, float CosAlpha, float SinAlpha, float x, float y, float z)
{
   float dum;
   float ax,ay,az;

   R[0]=R[5]=R[10]=R[15]=1.0;
   R[1]=R[2]=R[3]=0.0;
   R[4]=R[6]=R[7]=0.0;
   R[8]=R[9]=R[11]=0.0;
   R[12]=R[13]=R[14]=0.0;
      
   if(x==0.0 && y==0.0 && z==0.0) return;
   
   /* find the unit vector (ax,ay,az) in the direction of (x,y,z) */
   dum=(float)sqrt((double)x*x+y*y+z*z);
   ax=x/dum; ay=y/dum; az=z/dum;
      
   R[0] = ax*ax + CosAlpha - CosAlpha*ax*ax; 
   R[1] = ax*ay - CosAlpha*ax*ay - SinAlpha*az; 
   R[2] = ax*az - CosAlpha*ax*az + SinAlpha*ay; 
   
   R[4] = ay*ax - CosAlpha*ay*ax + SinAlpha*az; 
   R[5] = ay*ay + CosAlpha - CosAlpha*ay*ay; 
   R[6] = ay*az - CosAlpha*ay*az - SinAlpha*ax; 
   
   R[8]  = az*ax - CosAlpha*az*ax - SinAlpha*ay; 
   R[9]  = az*ay - CosAlpha*az*ay + SinAlpha*ax; 
   R[10] = az*az + CosAlpha - CosAlpha*az*az; 
   
   return;
}
