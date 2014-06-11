#include <babak_lib.h>

/* 
 This function determines the general direction of a vector (x,y,z) defined in RAS system.
 Possible directions are: Right, Left, Anterior, Posterior, Superior, or Inferior.
 Input: (x,y,z) a vector defined in RAS system
 Output: One of six charaters {R,L,A,P,S,I}
*/
char directionCode(float x, float y, float z) 
{
   float lxl, lyl, lzl; // absolute values of x, y, and z

   lxl = fabsf(x);
   lyl = fabsf(y);
   lzl = fabsf(z);

   if( lxl>lyl && lxl>lzl)
   {
      if(x>0.0) return('R'); else return('L');
   }
   else if( lyl>lxl && lyl>lzl)
   {
      if(y>0.0) return('A'); else return('P');
   }
   else 
   {
      if(z>0.0) return('S'); else return('I');
   }
}

