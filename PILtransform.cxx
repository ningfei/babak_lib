#include <babak_lib.h>

void PILtransform(const char *orientCode, float *orientMat)
{
   char c;

   for(int i=0; i<16; i++) orientMat[i]=0.0;

   orientMat[15]=1.0;

   for(int j=0; j<3; j++)
   {
      // read the jth element of the orientCode vector
      c = orientCode[j];		

      c=toupper(c);	// convert c to upper case

      // set the jth column of the orientMat
      if( c=='P' )		orientMat[j]   =  1.0;
      else if( c=='A' )	orientMat[j]   = -1.0;
      else if( c=='I' )	orientMat[4+j] =  1.0;
      else if( c=='S' )	orientMat[4+j] = -1.0;
      else if( c=='L' )	orientMat[8+j] =  1.0;
      else if( c=='R' )	orientMat[8+j] = -1.0;
   }

   return;
}

void inversePILtransform(const char *orientCode, float *orientMat)
{
   char c;

   for(int i=0; i<16; i++) orientMat[i]=0.0;

   orientMat[15]=1.0;

   for(int i=0; i<3; i++)
   {
      // read the ith element of the orientCode vector
      c = orientCode[i];		
   
      c=toupper(c);	// convert c to upper case

      // set the ith row of the orientMat
      if( c=='P' )		orientMat[4*i] 		=  1.0;
      else if( c=='A' )	orientMat[4*i] 		= -1.0;
      else if( c=='I' )	orientMat[4*i + 1]	=  1.0;
      else if( c=='S' )	orientMat[4*i + 1]	= -1.0;
      else if( c=='L' )	orientMat[4*i + 2]	=  1.0;
      else if( c=='R' )	orientMat[4*i + 2]	= -1.0;
   }

   return;
}
