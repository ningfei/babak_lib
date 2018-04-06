#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "babak_lib.h"

#define NO 0
#define YES 1
#define LIN 1

void print_help_and_exit()
{
	printf("Usage: multiplyTransformation <Tin1> <Tin2> <Tout>\n"
	"Computes Tout = Tin1*Tin2\n"
	);
	exit(0);
}

int main(int argc, char **argv)
{
	float Tin1[16],Tin2[16],Tout[16];

	printf("%d\n",argc);

	if(argc != 4)
		print_help_and_exit();

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// ensure the input transformatin file exists, is readable, and has the expected size
	if ( checkFileExistence(argv[1])==0 )
	{
		printf("Error: File %s does not exist! Aborting ...\n",argv[1]);
		exit(0);
	}

	if ( checkFileReadOK(argv[1])==0 )
	{
		printf("Error: Read permission for %s denied! Aborting ...\n",argv[1]);
		exit(0);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////

	loadTransformation( argv[1], Tin1);

	printf("\nTransformation matrix 1:\n");
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin1[0],Tin1[1],Tin1[2],Tin1[3]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin1[4],Tin1[5],Tin1[6],Tin1[7]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin1[8],Tin1[9],Tin1[10],Tin1[11]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin1[12],Tin1[13],Tin1[14],Tin1[15]);
	printf("\n");

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// ensure the input transformatin file exists, is readable, and has the expected size
	if ( checkFileExistence(argv[2])==0 )
	{
		printf("Error: File %s does not exist! Aborting ...\n",argv[2]);
		exit(0);
	}

	if ( checkFileReadOK(argv[2])==0 )
	{
		printf("Error: Read permission for %s denied! Aborting ...\n",argv[2]);
		exit(0);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////

	loadTransformation( argv[2], Tin2);

	printf("\nTransformation matrix 2:\n");
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin2[0],Tin2[1],Tin2[2],Tin2[3]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin2[4],Tin2[5],Tin2[6],Tin2[7]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin2[8],Tin2[9],Tin2[10],Tin2[11]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tin2[12],Tin2[13],Tin2[14],Tin2[15]);
	printf("\n");

	multi(Tin1,4,4,Tin2,4,4,Tout);

	printf("\nTout = Tin1 x Tin2\n");
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tout[0],Tout[1],Tout[2],Tout[3]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tout[4],Tout[5],Tout[6],Tout[7]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tout[8],Tout[9],Tout[10],Tout[11]);
	printf("%6.3f  %6.3f  %6.3f  %6.3f\n",Tout[12],Tout[13],Tout[14],Tout[15]);
	printf("\n");

	{
		FILE *fp;
		fp=fopen(argv[3],"w");
		for(int i=0; i<4; i++)
		{
				for(int j=0; j<4; j++)
				{
						fprintf(fp,"%f\t",Tout[i*4+j]);
				}
				fprintf(fp,"\n");
		}
		fclose(fp);
	}
}
