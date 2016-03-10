#include <stdlib.h>
#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "../include/spm_analyze.h"
#include "../include/babak_lib.h"
#include <f2c.h>


int main(int argc, char **argv)
{

/*
*/
	// this section tests the paired samples t-test routine 
	short x1[12];
	short x2[12];
	int n=7;
	int df;
	double meandiff,t;

	x1[0] = 78;
	x1[1] = 77;
	x1[2] = 70;
	x1[3] = 88;
	x1[4] = 85;
	x1[5] = 83;
	x1[6] = 80;

	printf("\nSample Mean 1 = %lf\n", sample_mean(x1,n));

	x2[0] = 83;
	x2[1] = 83;
	x2[2] = 85;
	x2[3] = 88;
	x2[4] = 89;
	x2[5] = 82;
	x2[6] = 77;

	printf("\nSample Mean 2 = %lf\n", sample_mean(x2,n));

	t=paired_samples_t(x1,x2,n,&df,&meandiff);
	printf("\nPaired samples t-value = %lf mean diff. = %lf df=%d\n",t,meandiff,df);

/*
	// this section tests the independent samples t-test routine 
	float x1[7];
	float x2[6];
	int n1=7;
	int n2=6;
	int df;
	double meandiff,t;

	x1[0] = 78;
	x1[1] = 77;
	x1[2] = 70;
	x1[3] = 88;
	x1[4] = 85;
	x1[5] = 83;
	x1[6] = 80;

	printf("\nSample Mean 1 = %lf\n", sample_mean(x1,n1));

	x2[0] = 83;
	x2[1] = 83;
	x2[2] = 85;
	x2[3] = 88;
	x2[4] = 89;
	x2[5] = 82;

	printf("\nSample Mean 2 = %lf\n", sample_mean(x2,n2));

	t=independent_samples_t(x1,n1,x2,n2,&df,&meandiff);
	printf("\nPaired samples t-value = %lf mean diff. = %lf\n",t,meandiff);
*/
}
