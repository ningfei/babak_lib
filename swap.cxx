#define _swap

#include <babak_lib.h>
#include <stdlib.h>
#include <stdio.h>

int bigEndian();
void swapN(char *in, int N);

////////////////////////////////////////////////////////////////////////////////
// Author: Babak A. Ardekani, Ph.D.
// Last update: April 12, 2010.
////////////////////////////////////////////////////////////////////////////////
void swapByteOrder(char *in, int N)
{
   char *dum;

   if(N<=0)
   {
      printf("\n\nWarning swapByteOrder(): Non-positive array dimension argument.\n\n");
      return;
   }

   if(in==NULL)
   {
      printf("\n\nWarning swapByteOrder(): NULL array argument.\n\n");
      return;
   }

   dum = (char *)malloc(N);

   if(dum==NULL)
   {
      printf("\n\nWarning swapByteOrder(): Memory allocation failure.\n\n");
      return;
   }

   for( int i=0; i<N; i++)
   {
      dum[i]=in[N-1-i];
   }

   for( int i=0; i<N; i++)
   {
      in[i]=dum[i];
   }

   free(dum);
}       
////////////////////////////////////////////////////////////////////////////////

// Returns 1 if the computer is bigEndian (e.g., SUN) and 0 if it is Little Endian (e.g., IBM PC)
int bigEndian()
{
	short s;
	char *cp;

	// We will set s=1. If the computer is Big Endian, the first byte in s will be
	// 0 and the second byte will be 1.  So this function returns the second byte.

	s=1;

	cp=(char *)(&s);
	return( (int)cp[1] );
}  

void swapN(char *in, int N)
{
	char dum[2]; 

	for(int i=0; i<N/2; i++)
	{    
		dum[0]=in[i*2+1];
		dum[1]=in[i*2];
        
		in[i*2]=dum[0];
		in[i*2+1]=dum[1];
	}
}       

void swap_float_array( float *x, int n)
{
   for(int i=0; i<n; i++)
   {    
      swapByteOrder( (char *)(x+i), sizeof(float) );
   }
}       

void swap_int_array( int *x, int n)
{
   for(int i=0; i<n; i++)
   {    
      swapByteOrder( (char *)(x+i), sizeof(int) );
   }
}       

void swap_model_file_hdr(model_file_hdr *hdr)
{
   swapByteOrder( (char *)(&hdr->nxHR), sizeof(int));
   swapByteOrder( (char *)(&hdr->nzHR), sizeof(int));
   swapByteOrder( (char *)(&hdr->dxHR), sizeof(float));
   swapByteOrder( (char *)(&hdr->nxLR), sizeof(int));
   swapByteOrder( (char *)(&hdr->nvol), sizeof(int));
   swapByteOrder( (char *)(&hdr->RPtemplateradius), sizeof(int));
   swapByteOrder( (char *)(&hdr->RPtemplateheight), sizeof(int));
   swapByteOrder( (char *)(&hdr->RPtemplatesize), sizeof(int));
   swapByteOrder( (char *)(&hdr->ACtemplateradius), sizeof(int));
   swapByteOrder( (char *)(&hdr->ACtemplateheight), sizeof(int));
   swapByteOrder( (char *)(&hdr->ACtemplatesize), sizeof(int));
   swapByteOrder( (char *)(&hdr->PCtemplateradius), sizeof(int));
   swapByteOrder( (char *)(&hdr->PCtemplateheight), sizeof(int));
   swapByteOrder( (char *)(&hdr->PCtemplatesize), sizeof(int));
   swapByteOrder( (char *)(&hdr->nangles), sizeof(int));
}

void swap_model_file_tail(model_file_tail *tail)
{
   swapByteOrder( (char *)(&tail->RPPCmean[0]), sizeof(float));
   swapByteOrder( (char *)(&tail->RPPCmean[1]), sizeof(float));
   swapByteOrder( (char *)(&tail->parcomMean), sizeof(float));
   swapByteOrder( (char *)(&tail->percomMean), sizeof(float));
   swapByteOrder( (char *)(&tail->RPmean[0]), sizeof(float));
   swapByteOrder( (char *)(&tail->RPmean[1]), sizeof(float));
}
