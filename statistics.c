#define _statistics

#include <stdlib.h> 
#include <math.h>
#include <stdio.h>
#include "babak_lib.h"

//////////////////////////////////////////////////////////////////
float imageMean(short *im, short *msk, int nv)
{
   int nbv=0;
   double sum=0.0;

   for(int i=0; i<nv; i++)
   {
      if(msk[i]>0)
      {
         nbv++;
         sum += im[i];
      }
   }

   if(nbv>0)
   {
      sum /= nbv;
   }

   return( (float)sum);
}

//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
float median(float *x, char *mask, int n)
{
   float *masked_x; // array x after removing the elements where mask[i]=0
   float m; // computed median
   unsigned long size; // size of masked_x array

   if(n<=0)
   {
      printf("\nWarning: invalid parameter n passed to median() function.\n");
      return(0.0);
   }

   if(n==1)
   {
      return(x[0]);
   }

   // since size<=n, we allocat an array of size n
   masked_x = (float *)calloc(n, sizeof(float));
   if( masked_x==NULL)
   {
      printf("\nError: Memory allocation failure for variable masked_x in median() function.\n");
      exit(1);
   }

   // if a NULL value is passed for mask, we just copy x to masked_x
   if(mask==NULL)
   {
      size = n;

      for(int i=0; i<n; i++)
      {
         masked_x[i] = x[i];
      }
   }
   else
   {
      size=0;
      for(int i=0; i<n; i++)
      {
         if( mask[i] )
         {
            masked_x[size] = x[i];
            size++;
         }
      }
   }

   // if all elements of mask are 0, pretend there is no mask
   if(size==0) 
   {
      printf("\nWarning: invalid mask passed to median() function.\n");
      size = n;

      for(int i=0; i<n; i++)
      {
         masked_x[i] = x[i];
      }
   }

   hpsort(size, masked_x);

   if( size%2 )
   {
      m = masked_x[ (size-1)/2 ];
   }
   else
   {
      m = 0.5*(masked_x[ size/2 ] + masked_x[ size/2 - 1]);
   }

   free(masked_x);

   return(m);
}
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
void decomposeVector(double *x, double *xpar, double *xper, double *u, int n)
{
	double c;

	normalize(u,n);

	c = dot(x,u,n);

	for(int i=0; i<n; i++) 
	{
		xpar[i] = c*u[i];
		xper[i] = x[i] - c*u[i];
	}
}


//////////////////////////////////////////////////////////////////
// Removes the component of vector x that is parrallel to vector u.
// Vector u will be normalized (i.e., ||u|| = 1.0)
void removeComponent(double *x, double *u, int n)
{
	double c;

	normalize(u,n);

	c = dot(x,u,n);

	for(int i=0; i<n; i++) x[i] -= c*u[i];
}

void removeComponent(float *x, float *u, int n)
{
	float c;

	normalize(u,n);

	c = dot(x,u,n);

	for(int i=0; i<n; i++) x[i] -= c*u[i];
}

/////////////////////////////////////////////////////////
// y is a (nxp) matrix.  This function operates on the columns of y.
// Upon completion, each of the p columns of y will be scale such that
// the maximum absolute value in each column is 1.
void scaleAbsToOne(double *y, int n, int p)
{
	double val;

	if(n<=0 || p<=0) return;

	for(int j=0; j<p; j++)
	{
		val = 0.0;

		for(int i=0; i<n; i++)
		if( fabs(y[i*p +j]) > val ) val = fabs(y[i*p + j]);

		if(val != 0.0)
		for(int i=0; i<n; i++)
			y[i*p +j] = y[i*p +j ]/val;
	}
}

void scaleAbsToOne(float *y, int n, int p)
{
	float val;

	if(n<=0 || p<=0) return;

	for(int j=0; j<p; j++)
	{
		val = 0.0;

		for(int i=0; i<n; i++)
		if( fabsf(y[i*p +j]) > val ) val = fabsf(y[i*p + j]);

		if(val != 0.0)
		for(int i=0; i<n; i++)
			y[i*p +j] = y[i*p +j ]/val;
	}
}

// y is an n-vector.  This function operates on y.
// Upon completion, y will be scaled such that
// the maximum absolute value is 1.
double scaleAbsToOne(double *y, int n)
{
   double val;

   if(n<=0) return(1.0);

   val = 0.0;

   for(int i=0; i<n; i++)
      if( fabs(y[i]) > val ) val = fabs(y[i]);

   if(val != 0.0)
   for(int i=0; i<n; i++)
      y[i] = y[i]/val;

   return(val);
}
/////////////////////////////////////////////////////////

// y is a (nxp) matrix.  This function operates on the columns of y.
// Upon completion, each of the p columns of y will have zero mean.
void removeVectorMean(double *y, int n, int p)
{
	double mean;

	if(n<=0 || p<=0) return;

	for(int j=0; j<p; j++)
	{
		mean = 0.0;

		for(int i=0; i<n; i++)
			mean += y[i*p + j];

		mean /= n;

		for(int i=0; i<n; i++)
			y[i*p +j] = y[i*p +j ] - mean;
	}
}

// y is a (nxp) matrix.  This function operates on the columns of y.
// Upon completion, each of the p columns of y will have zero mean.
void removeVectorMean(float *y, int n, int p)
{
	float mean;

	if(n<=0 || p<=0) return;

	for(int j=0; j<p; j++)
	{
		mean = 0.0;

		for(int i=0; i<n; i++)
			mean += y[i*p + j];

		mean /= n;

		for(int i=0; i<n; i++)
			y[i*p +j] = y[i*p +j ] - mean;
	}
}

double removeVectorMean(double *x, double *y, int n)
{
	double mean=0.0;

	if(n<=0) return(0.0);

	for(int i=0; i<n; i++)
		mean += x[i];

	mean /= n;

	for(int i=0; i<n; i++)
		y[i] = x[i]-mean;

	return(mean);
}

double removeVectorMean(short *x, double *y, int n)
{
	double mean=0.0;

	if(n<=0) return(0.0);

	for(int i=0; i<n; i++)
		mean += x[i];

	mean /= n;

	for(int i=0; i<n; i++)
		y[i] = x[i]-mean;

	return(mean);
}

double removeVectorMean(double *y, int n)
{
	double mean=0.0;

	if(n<=0) return(0.0);

	for(int i=0; i<n; i++)
		mean += y[i];

	mean /= n;

	for(int i=0; i<n; i++)
		y[i] = y[i]-mean;

	return(mean);
}

double removeVectorMean(float *x, float *y, int n)
{
	double mean=0.0;

	if(n<=0) return(0.0);

	for(int i=0; i<n; i++)
		mean += x[i];

	mean /= n;

	for(int i=0; i<n; i++)
		y[i] = (float)(x[i]-mean);

	return(mean);
}

double removeVectorMean(float *y, int n)
{
	double mean=0.0;

	if(n<=0) return(0.0);

	for(int i=0; i<n; i++)
		mean += y[i];

	mean /= n;

	for(int i=0; i<n; i++)
		y[i] = (float)(y[i]-mean);

	return(mean);
}

/////////////////////////////////////////////////////////

// Implements Eq. (3.3.11) of J. Cohen & P. Cohen (2nd ed.)
void partialCorrelation(float *Y, float *X1, float *X2, int n, double *pr1, double *pr2)
{
	double rY1, rY2, r12;
	double dum1, dum2;

	rY1 = pearsonCorrelation(X1,Y, n);
	rY2 = pearsonCorrelation(X2,Y, n);
	r12 = pearsonCorrelation(X1, X2, n);

	dum1 = (1.0-rY2*rY2)*(1-r12*r12);
	dum2 = (1.0-rY1*rY1)*(1-r12*r12);

	if(dum1>0.0) 
		*pr1 = (rY1 - rY2*r12)/sqrt(dum1);
	else 
		*pr1 = 0.0;

	if(dum2>0.0) 
		*pr2 = (rY2 - rY1*r12)/sqrt(dum2);
	else
		*pr2 = 0.0;
}

void partialCorrelation(double *Y, double *X1, double *X2, int n, double *pr1, double *pr2)
{
	double rY1, rY2, r12;
	double dum1, dum2;

	rY1 = pearsonCorrelation(X1,Y, n);
	rY2 = pearsonCorrelation(X2,Y, n);
	r12 = pearsonCorrelation(X1, X2, n);

	dum1 = (1.0-rY2*rY2)*(1-r12*r12);
	dum2 = (1.0-rY1*rY1)*(1-r12*r12);

	if(dum1>0.0) 
		*pr1 = (rY1 - rY2*r12)/sqrt(dum1);
	else 
		*pr1 = 0.0;

	if(dum2>0.0) 
		*pr2 = (rY2 - rY1*r12)/sqrt(dum2);
	else
		*pr2 = 0.0;
}

// Implements Eq. (3.3.11) of J. Cohen & P. Cohen (2nd ed.)
void partialCorrelation(short *Y, float *X1, float *X2, int n, double *pr1, double *pr2)
{
	double rY1, rY2, r12;
	double dum1, dum2;

	rY1 = pearsonCorrelation(X1, Y, n);
	rY2 = pearsonCorrelation(X2, Y, n);
	r12 = pearsonCorrelation(X1, X2, n);

	dum1 = (1.0-rY2*rY2)*(1-r12*r12);
	dum2 = (1.0-rY1*rY1)*(1-r12*r12);

	if(dum1>0.0) 
		*pr1 = (rY1 - rY2*r12)/sqrt(dum1);
	else 
		*pr1 = 0.0;

	if(dum2>0.0) 
		*pr2 = (rY2 - rY1*r12)/sqrt(dum2);
	else
		*pr2 = 0.0;
}

////////////////////////////////////////////////////////////////

double one_sample_t(double *x, int n)
{
   double var;
   double mean;
   double t;

   if(n<2) 
   {
      return(0.0);
   }

   var = sample_variance(x,n,mean);

   if(var<=0.0) 
   {
      return(0.0);
   }

   t = sqrt(1.0*n)*mean/sqrt(var);

   return(t);
}
