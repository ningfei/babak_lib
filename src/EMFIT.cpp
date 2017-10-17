#define _EMFIT

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void EMFIT1d(double *hist, double *fit, short *label, int nb, double *mean, double *var, double *p, int nclass, int niter)
{
   double pi;
   double val;
   double *g;
   double sum=0.0;

   pi = 4.0*atan(1.0);

   g = (double *)calloc(nclass*nb,sizeof(double));

   // initialization
   ///////////////////////////////////////////////////////////////
   // assigns bins to different classes
   for(int b=0; b<nb; b++)
   {
      label[b] = (b*nclass)/nb;
   }

   for(int k=0; k<nclass; k++)  // for each row
   {
      for(int b=0; b<nb; b++) 
      {
         if( label[b]==k )
            g[k*nb + b]=hist[b];
         else
            g[k*nb + b]=0.0;
      }
   }
   ///////////////////////////////////////////////////////////////

   for(int iter=0; iter<niter; iter++)
   {
      for(int k=0; k<nclass; k++) 
      {
         p[k] = 0.0;

         for(int b=0; b<nb; b++) p[k] += g[k*nb + b];
      }

      sum = 0.0;
      for(int k=0; k<nclass; k++) sum += p[k];

      if( sum > 0.0) for(int k=0; k<nclass; k++) p[k]/=sum;

      for(int k=0; k<nclass; k++) if(p[k]<0.0001) p[k]=0.0;

      for(int k=0; k<nclass; k++) 
      if(p[k]>0.0)
      {
         mean[k] = 0.0;

         for(int b=0; b<nb; b++) mean[k] += b*g[k*nb + b];

         mean[k] /= p[k];
      }
      else mean[k] = 0.0;

		for(int k=0; k<nclass; k++) 
		if(p[k]>0.0)
		{
			var[k] = 0.0;

			for(int b=0; b<nb; b++) var[k] += (b-mean[k])*(b-mean[k])*g[k*nb + b];

			var[k] /= p[k];
		}
		else	 var[k] = 0.0;

		for(int b=0; b<nb; b++) 
		{
			double sum = 0.0;

			for(int k=0; k<nclass; k++) 
			if(var[k]>0.0)
			{
				val = p[k]*exp(-0.5*(b-mean[k])*(b-mean[k])/var[k])/sqrt(var[k]*2*pi);
				g[k*nb+b] = hist[b]*val;
				sum += val;
			}
			else	 g[k*nb+b] = 0.0;

			for(int k=0; k<nclass; k++) 
			if(sum>0.0) g[k*nb+b] /= sum;
		}
	}

	for(int b=0; b<nb; b++) 
	{
		double maxval = 0.0;
		for(int k=0; k<nclass; k++) 
		if(p[k]>0.0)
		{
			if(g[k*nb+b]>maxval) { maxval=g[k*nb+b]; label[b]=k; }
		}
	}

	for(int i=0; i<nb; i++) fit[i]=0.0;

	for(int b=0; b<nb; b++)
	{
		for(int k=0; k<nclass; k++) 
		if(p[k]>0.0 && var[k]>0.0)
		{
			fit[b] += exp(-0.5*(b-mean[k])*(b-mean[k])/var[k])*p[k]/sqrt(var[k]*2*pi);
		}
	}

	free(g);
}
