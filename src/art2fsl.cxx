#include <stdlib.h>
#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include "../include/spm_analyze.h"
#include "../include/babak_lib.h"
#include "../include/sph.h"
#include "../include/landmarks.h"
#include "../include/minmax.h"
#include <ctype.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
   {"-i", 1, 'i'},
   {"-o", 1, 'o'},
   {"-v", 0, 'v'},
   {"-h", 0, 'h'},
   {"-Tdx",1,'X'},
   {"-Tdy",1,'Y'},
   {"-Tdz",1,'Z'},
   {"-Tnx",1,'x'},
   {"-Tny",1,'y'},
   {"-Tnz",1,'z'},
   {"-Sdx",1,'1'},
   {"-Sdy",1,'2'},
   {"-Sdz",1,'3'},
   {"-Snx",1,'4'},
   {"-Sny",1,'5'},
   {"-Snz",1,'6'},
   {"-help", 0, 'h'},
   {0, 0, 0}
};

void print_help_and_exit()
{
   printf("\nUsage: art2fsl [-v -h] -i <ART matrix> -o <FSL matrix>\n"
   "-Tnx <n> -Tny <n> -Tnz <n> -Tdx <f> -Tdy <f> -Tdz <f>\n"
   "-Snx <n> -Sny <n> -Snz <n> -Sdx <f> -Sdy <f> -Sdz <f>\n"
   "Required:\n"
   "\t-i <image list file>: A text file containing the list of NIFTI images to be registered\n"
   "\t-nx <int> -ny <int> -nz <int>: Output matrix size (default: 255x255x189)\n\n"
   "\t-dx <float> -dy <float> -dz <float>: Output voxel size (default: 1.0x1.0x1.0 mm^3)\n\n"
   "Options:\n"
   "\t-v Enables verbose mode\n\n" 
   );

   exit(0);
}

int main(int argc, char **argv)
{
   char inputmatrixfile[DEFAULT_STRING_LENGTH]="";
   char outputmatrixfile[DEFAULT_STRING_LENGTH]="";
   float Mart[16];
   float Mfsl[16];
   DIM sub_dim, trg_dim;
   FILE *fp;

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case '1':
            sub_dim.dx = atof(optarg);
            break;
         case '2':
            sub_dim.dy = atof(optarg);
            break;
         case '3':
            sub_dim.dz = atof(optarg);
            break;
         case '4':
            sub_dim.nx = atoi(optarg);
            break;
         case '5':
            sub_dim.ny = atoi(optarg);
            break;
         case '6':
            sub_dim.nz = atoi(optarg);
            break;
         case 'X':
            trg_dim.dx = atof(optarg);
            break;
         case 'Y':
            trg_dim.dy = atof(optarg);
            break;
         case 'Z':
            trg_dim.dz = atof(optarg);
            break;
         case 'x':
            trg_dim.nx = atoi(optarg);
            break;
         case 'y':
            trg_dim.ny = atoi(optarg);
            break;
         case 'z':
            trg_dim.nz = atoi(optarg);
            break;
         case 'i':
            sprintf(inputmatrixfile,"%s",optarg);
            break;
         case 'o':
            sprintf(outputmatrixfile,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '?':
            print_help_and_exit();
      }
   }

   if(opt_v)
   {
      printf("Input (ART) transformation matrix file: %s\n",inputmatrixfile);
      printf("Output (FSL) transformation matrix file: %s\n",outputmatrixfile);
      printf("Target (reference) image matrix dimensions: %d x %d x %d\n",trg_dim.nx,trg_dim.ny,trg_dim.nz);
      printf("Target (reference) image voxel dimensions: %6.4f x %6.4f x %6.4f\n",trg_dim.dx,trg_dim.dy,trg_dim.dz);
      printf("Subject (moving) image matrix dimensions: %d x %d x %d\n",sub_dim.nx,sub_dim.ny,sub_dim.nz);
      printf("Subject (moving) image voxel dimensions: %6.4f x %6.4f x %6.4f\n",sub_dim.dx,sub_dim.dy,sub_dim.dz);
   }

   loadTransformation(inputmatrixfile, Mart);

   if(opt_v)
   {
      printMatrix(Mart, 4, 4, "Input ART transformation:", NULL);
   }

   art_to_fsl(Mart, Mfsl, sub_dim, trg_dim);

   if(opt_v)
   {
      printMatrix(Mfsl, 4, 4, "Output FSL transformation:", NULL);
   }

   fp = fopen(outputmatrixfile,"w");
   fprintf(fp,"%f %f %f %f\n",Mfsl[0],Mfsl[1],Mfsl[2],Mfsl[3]);
   fprintf(fp,"%f %f %f %f\n",Mfsl[4],Mfsl[5],Mfsl[6],Mfsl[7]);
   fprintf(fp,"%f %f %f %f\n",Mfsl[8],Mfsl[9],Mfsl[10],Mfsl[11]);
   fprintf(fp,"%f %f %f %f\n",Mfsl[12],Mfsl[13],Mfsl[14],Mfsl[15]);
   fclose(fp);

   return 0;
}

