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
#include "spm_analyze.h"
#include "babak_lib.h"
#include <f2c.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
        {"-cc", 1, 'c'},
        {"-t1", 1, '1'},
        {"-t2", 1, '2'},
        {"-o", 1, 'o'},
        {"-i", 1, 'i'},
        {0, 0, 0}
};

int opt_log=NO;
int opt_t=NO;
int opt_d=NO;

int opt_o=NO;
int opt_p=NO;
int opt_P=NO;
int opt_cc=NO;

void print_help_and_exit()
{
	exit(0);
}

int main(int argc, char **argv)
{
   int ccthresh=0;

   nifti_1_header hdr;

	int ncc, ntcc;

	char inputfile[512];
	char outputfile[1024];
	float *ip_image;
	float thresh1, thresh2;

	char dfFile[512];
	short *dfmap;
	double df;

	int nx, ny, nz, np, nv;
	float dx, dy, dz;

	FILE *logFilePtr;

	int count=0;
	char logFile[1024];

	char prefix[1024];
	double alpha, alpha2;

	while( (opt=getoption(argc, argv, options)) != -1)
	{
		switch (opt) {
			case 'c':
				ccthresh=atoi(optarg);
				break;
			case 'i':
				sprintf(inputfile,"%s",optarg);
				break;
			case 'o':
				sprintf(outputfile,"%s",optarg);
				break;
			case '1':
				thresh1 = atof(optarg);
				break;
			case '2':
				thresh2 = atof(optarg);
				break;
			case '?':
				print_help_and_exit();
		}
	}

   if(argc==1) print_help_and_exit();

   printf("CC threshold = %d\n",ccthresh);
   printf("Thresh 1 = %f\n",thresh1);
   printf("Thresh 2 = %f\n",thresh2);

   ip_image = (float *)read_nifti_image(inputfile, &hdr);
   nx = hdr.dim[1]; ny = hdr.dim[2]; nz = hdr.dim[3];
   dx = hdr.pixdim[1]; dy = hdr.pixdim[2]; dz = hdr.pixdim[3];

   nv = nx*ny*nz;

	for(int i=0; i<nv; i++)
    {
	   if(ip_image[i]>=0 && ip_image[i]<=thresh1)
         ip_image[i]=0.0;
	   else if(ip_image[i]<0 && ip_image[i]>=thresh2)
         ip_image[i]=0.0;
    }

   if(ccthresh>0) 
   {
      char *im;
      int ncc, ntcc;

      im = (char *)calloc(nv, sizeof(char));

      for(int i=0; i<nv; i++) if(ip_image[i]!=0.0) im[i]=1; else im[i]=0;

      thr_Connected_Component(im, ccthresh, nx, ny, nz, &ncc, &ntcc);

      for(int i=0; i<nv; i++) if(im[i]==0) ip_image[i]=0.0;

      free(im);

      printf("\n%d clusters of %d had size >= %d\n",ntcc,ncc,ccthresh);
   }

   save_nifti_image(outputfile, ip_image, &hdr);
}
