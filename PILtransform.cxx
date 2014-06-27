#include <babak_lib.h>

void PILtransform(const char *inputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<16; i++) T[i]=0.0;
   T[15]=1.0;

   for(int j=0; j<3; j++)
   {
      // read the jth element of the inputOrientCode vector
      c = inputOrientCode[j];		

      c=toupper(c);	// convert c to upper case

      // set the jth column of the T
      if( c=='P' )		T[j]   =  1.0;
      else if( c=='A' )	T[j]   = -1.0;
      else if( c=='I' )	T[4+j] =  1.0;
      else if( c=='S' )	T[4+j] = -1.0;
      else if( c=='L' )	T[8+j] =  1.0;
      else if( c=='R' )	T[8+j] = -1.0;
   }

   return;
}

void inversePILtransform(const char *outputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<16; i++) T[i]=0.0;
   T[15]=1.0;

   for(int i=0; i<3; i++)
   {
      // read the ith element of the outputOrientCode vector
      c = outputOrientCode[i];		
   
      c=toupper(c);	// convert c to upper case

      // set the ith row of the T
      if( c=='P' )		T[4*i] 		=  1.0;
      else if( c=='A' )	T[4*i] 		= -1.0;
      else if( c=='I' )	T[4*i + 1]	=  1.0;
      else if( c=='S' )	T[4*i + 1]	= -1.0;
      else if( c=='L' )	T[4*i + 2]	=  1.0;
      else if( c=='R' )	T[4*i + 2]	= -1.0;
   }

   return;
}
