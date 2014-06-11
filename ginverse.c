#define _ginverse

#include <stdlib.h>
#include <stdio.h>
#include "babak_lib.h"


// X is an Nxp matrix
// G is a pxp matrix
// This function computes a 'generalized inverse' of matrix X'X and returns it in G.
// The generalized inverse is such that X'X * G * X'X = X'X. The return value of
// this function is the rank of matrix X'X. The return value is -1 if error occurs.

int ginverse(float *X, int N, int p, float *G)
{
	float *XtX;		// X'X
	float *D;
	float *Ut;
	float *DUt;
	int r;			// rank of XtX

	// ensure memory was allocated for X and G
	if(X==NULL || G==NULL)
	{
		printf("\n\nginverse(): NULL input matrix encountered.\n\n");
		return(-1);
	}

	// make sure N and p are positive
	if(N<=0 || p<=0)
	{
		printf("\n\nginverse(): Non-positive matrix dimension encountered.\n\n");
		return(-1);
	}

	// allocate memory of the matrix X'X
	XtX=(float *)calloc(p*p,sizeof(float));
	if(XtX==NULL) 
	{
		printf("\n\nginverse(): Memory allocation problem.\n\n");
		return(-1);
	}

	// compute X'X
	// Only the lower part of X'X is formed since it is symmetric.
	for(int i=0; i<p; i++)
	for(int j=0; j<=i; j++)
	{
		XtX[i*p+j]=0.0;

		// Element (i,j) of X'X is the dot product of columns i and j of X
		for(int k=0; k<N; k++) XtX[i*p+j] += X[k*p + i]*X[k*p + j];	
	}

	// Diagnolized X'X,	such that X'X = U diag(D) U'
	// This function returns U' in place of X'X
	D=diagATA_float(XtX, p, 'L');
	Ut=XtX;

	if(D==NULL)
	{
		printf("\n\nginverse(): Could not diagnolize X'X.\n\n");
		return(-1);
	}

	// Determine rank of X'X based on it's number of positive eigenvalues
	r=0;
	for(int i=0; i<p; i++) 
	{
		if( (int)(1000.0*D[i]/D[0])>0) 
		//if( D[i]>0.0 ) 
			r++; 
		else 
			D[i]=0.0;
	}

//	printf("%f %f %f\n",D[0],D[1],D[2]);

	// Set the last (p-r) rows of Ut to 0.
	for(int i=r; i<p; i++)
	for(int j=0; j<p; j++) Ut[i*p + j]=0.0;

	DUt=(float *)calloc(p*p,sizeof(float));
	if(DUt==NULL) 
	{
		printf("\n\nginverse(): Memory allocation problem.\n\n");
		return(-1);
	}

	// compute diag(1/D) * U'
	for(int i=0; i<r; i++)
	for(int j=0; j<p; j++) DUt[i*p + j]=Ut[i*p +j]/D[i];

	// Compute G = U * diag(1/D) * U'
	mat_trans_mat(Ut, p, p, DUt, p, G);

	free(XtX);
	free(DUt);
	free(D);

	return(r);
}

// X is an Nxp matrix
// G is a pxp matrix
// This function computes a 'generalized inverse' of matrix X'X and returns it in G.
// The generalized inverse is such that X'X * G * X'X = X'X. The return value of
// this function is the rank of matrix X'X. The return value is -1 if error occurs.

int ginverse(double *X, int N, int p, float *G)
{
	float *XtX;		// X'X
	float *D;
	float *Ut;
	float *DUt;
	int r;			// rank of XtX

	// ensure memory was allocated for X and G
	if(X==NULL || G==NULL)
	{
		printf("\n\nginverse(): NULL input matrix encountered.\n\n");
		return(-1);
	}

	// make sure N and p are positive
	if(N<=0 || p<=0)
	{
		printf("\n\nginverse(): Non-positive matrix dimension encountered.\n\n");
		return(-1);
	}

	// allocate memory of the matrix X'X
	XtX=(float *)calloc(p*p,sizeof(float));
	if(XtX==NULL) 
	{
		printf("\n\nginverse(): Memory allocation problem.\n\n");
		return(-1);
	}

	// compute X'X
	// Only the lower part of X'X is formed since it is symmetric.
	for(int i=0; i<p; i++)
	for(int j=0; j<=i; j++)
	{
		XtX[i*p+j]=0.0;

		// Element (i,j) of X'X is the dot product of columns i and j of X
		for(int k=0; k<N; k++) XtX[i*p+j] += (float)(X[k*p + i]*X[k*p + j]);	
	}

	// Diagnolized X'X,	such that X'X = U diag(D) U'
	// This function returns U' in place of X'X
	D=diagATA_float(XtX, p, 'L');
	Ut=XtX;

	if(D==NULL)
	{
		printf("\n\nginverse(): Could not diagnolize X'X.\n\n");
		return(-1);
	}

	// Determine rank of X'X based on it's number of positive eigenvalues
	r=0;
	for(int i=0; i<p; i++) 
	{
		if( (int)(1000.0*D[i]/D[0])>0) 
		//if( D[i]>0.0 ) 
			r++; 
		else 
			D[i]=0.0;
	}

//	printf("%f %f %f\n",D[0],D[1],D[2]);

	// Set the last (p-r) rows of Ut to 0.
	for(int i=r; i<p; i++)
	for(int j=0; j<p; j++) Ut[i*p + j]=0.0;

	DUt=(float *)calloc(p*p,sizeof(float));
	if(DUt==NULL) 
	{
		printf("\n\nginverse(): Memory allocation problem.\n\n");
		return(-1);
	}

	// compute diag(1/D) * U'
	for(int i=0; i<r; i++)
	for(int j=0; j<p; j++) DUt[i*p + j]=Ut[i*p +j]/D[i];

	// Compute G = U * diag(1/D) * U'
	mat_trans_mat(Ut, p, p, DUt, p, G);

	free(XtX);
	free(DUt);
	free(D);

	return(r);
}
