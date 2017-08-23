//////////////////////////////////////////////
// Examle VSPS detection failure
// 141_S_4438_??.nii
//
// Example of failed AC detection
// 053_S_2396_B0.nii
//
// Examles of failed MSP detection
// 082_S_4208_F1.nii
// 
// Example of bad image
// 153_S_4125_F0.nii
//
// AC/PC/RP detection failure (interesting case)
// Shows robustness of 8-landmark detection
// 941_S_5124_*.nii
//
// Excellent example of registration showing
// the improvement over MSP/AC-PC/LM PIL transformation
//031_S_2022_B0.nii
//031_S_2022_B1.nii
//031_S_2022_F0.nii
//031_S_2022_F1.nii
//
// Good for making the case for removing non-brain regions
//032_S_4429_B0.nii
//032_S_4429_B1.nii
//032_S_4429_F0.nii
//032_S_4429_F1.nii
//
// Good example: 
// 011_S_4912_B0.nii
// 011_S_4912_B1.nii
// 011_S_4912_F0.nii
// 011_S_4912_F1.nii
//
// Excellent example:
// 002_S_4171_B0.nii
// 002_S_4171_B1.nii
// 002_S_4171_F0.nii
// 002_S_4171_F1.nii
//
// Example of robustness wrt motion
// 011_S_4845_B0.nii
// 011_S_4845_B1.nii
// 011_S_4845_F0.nii
// 011_S_4845_F1.nii
//
// Good example, also handles a little motion
// 012_S_4545_B0.nii
// 012_S_4545_B1.nii
// 012_S_4545_F0.nii
// 012_S_4545_F1.nii
//
// Example of the usefullness of the braincloud mask
// 014_S_4039_B0.nii
// 014_S_4039_B1.nii
// 014_S_4039_F0.nii
// 014_S_4039_F1.nii
//
// Good examle of why non-brain regions need to be masked
// 024_S_4223_B0.nii
// 024_S_4223_B1.nii
// 024_S_4223_F0.nii
// 024_S_4223_F1.nii
//
//////////////////////////////////////////////

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

#define YES 1
#define NO 0

#define SNX 255
#define SNY 255
#define SNZ 189
#define SDX 1.0
#define SDY 1.0
#define SDZ 1.0

#define MAXIM 256 // maximum number of images allowed 

#define DEFAULT_SEARCH_RADIUS 5
#define DEFAULT_PATCh_RADIUS 7
#define MAX_RADIUS 100
#define MAXITER 25  // default maxiter
#define DEFAULT_STRING_LENGTH 512

int delx=20;
int dely=20;
int delz=20;

int maxiter=MAXITER;

int opt;

static struct option options[] =
{
        {"-v", 0, 'v'},
        {"-imlist", 1, 'i'},
        {"-i", 1, 'i'},
        {"-maxiter", 1, 'm'},
        {"-threshold", 1, 't'},
        {"-r", 1, 'r'},
        {"-R", 1, 'R'},
        {"-o", 1, 'o'},
        {"-h", 0, 'h'},
        {"-help", 0, 'h'},
        {0, 0, 0}
};

int opt_o=NO;
int opt_b=NO;

void print_help_and_exit()
{
   printf("\nUsage: atra [-h -help -v -threshold <threshold>] -i <image list file>...>\n"
   "-v Enables verbose mode\n" 
   "-h or -help: Prints help message\n" 
   "-o <output image name>: Specifies the filename for the outputted average image\n\n"
   "-threshold <threshold>: Binarizes the average images using <threhsold> level\n\n"
   "\n\n");

   exit(0);
}

// returns 1 if all images have the same dimensions nx, ny, and nz, 0 otherwise
int checkDimension_avgImage(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz)
{
   nifti_1_header hdr;
   short dataType;

   if(N==0) return(1);

   printf("Image %d: %s\n",1,imagefile[0]);

   hdr = read_NIFTI_hdr(imagefile[0]);
   *nx = hdr.dim[1];
   *ny = hdr.dim[2];
   *nz = hdr.dim[3];
   *dx = hdr.pixdim[1];
   *dy = hdr.pixdim[2];
   *dz = hdr.pixdim[3];
   dataType = hdr.datatype;

   for(int i=1; i<N; i++)
   {
      printf("Image %d: %s\n",i+1,imagefile[i]);
      hdr = read_NIFTI_hdr(imagefile[i]);

      if( *nx != hdr.dim[1] ||  *ny != hdr.dim[2] ||  *nz != hdr.dim[3]) 
      {
            return(0);
      }
   }

   return(1);
}

void set_hdr(nifti_1_header &hdr)
{
      float T[16];
      float T_ijk2xyz[16];
      mat44 R;

      hdr.dim[1]=SNX;
      hdr.dim[2]=SNY;
      hdr.dim[3]=SNZ;
      hdr.pixdim[1]=SDX;
      hdr.pixdim[2]=SDY;
      hdr.pixdim[3]=SDZ;

      // T will take points from PIL to RAS space
      inversePILtransform("RAS", T);

      ijk2xyz(T_ijk2xyz, SNX, SNY, SNZ, SDX, SDY, SDZ);

      // Takes points from (i,j,k) to RAS
      multi(T, 4, 4,  T_ijk2xyz, 4,  4, T);

      hdr.sform_code = 3;
      hdr.srow_x[0]=T[0]; hdr.srow_x[1]=T[1]; hdr.srow_x[2]=T[2]; hdr.srow_x[3]=T[3];
      hdr.srow_y[0]=T[4]; hdr.srow_y[1]=T[5]; hdr.srow_y[2]=T[6]; hdr.srow_y[3]=T[7];
      hdr.srow_z[0]=T[8]; hdr.srow_z[1]=T[9]; hdr.srow_z[2]=T[10]; hdr.srow_z[3]=T[11];

      hdr.qform_code = 3;
      hdr.qform_code = 3;
      R.m[0][0]=T[0];  R.m[0][1]=T[1];  R.m[0][2]=T[2];  R.m[0][3]=T[3];
      R.m[1][0]=T[4];  R.m[1][1]=T[5];  R.m[1][2]=T[6];  R.m[1][3]=T[7];
      R.m[2][0]=T[8];  R.m[2][1]=T[9];  R.m[2][2]=T[10]; R.m[2][3]=T[11];
      R.m[3][0]=T[12]; R.m[3][1]=T[13]; R.m[3][2]=T[14]; R.m[3][3]=T[15];

      nifti_mat44_to_quatern( R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
      &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
      &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
}

void compute_landmark_cm(int nim, int lm_in[][3], int lmcm[], char loomsk[])
{
   float x=0.0, y=0.0, z=0.0;
   int n=0;

   for(int i=0; i<nim; i++)
   if( loomsk[i]==1 )
   {
      x += lm_in[i][0];
      y += lm_in[i][1];
      z += lm_in[i][2];
      n++;
   }

   if( n>0 )
   {
      lmcm[0] = (int)(x/n + 0.5);
      lmcm[1] = (int)(y/n + 0.5);
      lmcm[2] = (int)(z/n + 0.5);
   }
}

// The average of rows of A is returned in sph.v
void extract_patches(int nim, SPH &sph, SHORTIM im[], int lm[][3], float *A)
{
   for(int i=0; i<nim; i++)
   {
      sph.set(im[i], lm[i][0], lm[i][1], lm[i][2]); 
      sph.zeromean();
      sph.normalize();
      sph.get(A+i*sph.n);
   }

   avgRow(A, nim, sph.n, sph.v);
}

// The average of rows of A is returned in sph.v
void extract_patches(int nim, SPH &sph, SHORTIM im[], int lm_in[][3], float *A, char loomsk[])
{
   int nrow=0;

   for(int i=0; i<nim; i++)
   if( loomsk[i]==1 )
   {
      sph.set(im[i], lm_in[i][0], lm_in[i][1], lm_in[i][2]); 
      sph.zeromean();
      sph.normalize();
      sph.get(A + nrow*sph.n);
      nrow++;
   }

   avgRow(A, nrow, sph.n, sph.v);
}

int lm_distance(int lm1[][3], int lm2[][3], int n)
{
   int d=0;

   for(int i=0; i<n; i++)
   {
      if( lm1[i][0] > lm2[i][0] ) d += (lm1[i][0]-lm2[i][0]);
      else d += (lm2[i][0]-lm1[i][0]);

      if( lm1[i][1] > lm2[i][1] ) d += (lm1[i][1]-lm2[i][1]);
      else d += (lm2[i][1]-lm1[i][1]);

      if( lm1[i][2] > lm2[i][2] ) d += (lm1[i][2]-lm2[i][2]);
      else d += (lm2[i][2]-lm1[i][2]);
   }

   return(d);
}

void lm_copy(int lm1[][3], int lm2[][3], int n)
{
   for(int i=0; i<n; i++)
   {
      lm2[i][0]=lm1[i][0];
      lm2[i][1]=lm1[i][1];
      lm2[i][2]=lm1[i][2];
   }
}

// nim: number of images in the training set
// refsph: spherical patch holder
int seek_lm(int nim, SPH &refsph, SHORTIM im[], int lm_in[][3], float *A, float *ssd, char loomsk[], int lmcm[],
SPH &searchsph, SPH &testsph, int lm_out[][3])
{
   int d;
   int search_center[3];

   search_center[0]=lm_in[0][0];
   search_center[1]=lm_in[0][1];
   search_center[2]=lm_in[0][2];

   for(int i=0; i<nim; i++) loomsk[i]=1;

   for(int iter=0; iter<maxiter; iter++)
   {
      float dum=0.0; 
      extract_patches(nim, refsph, im, lm_in, A);
      ssdRow(A, nim, refsph.n, refsph.v, ssd);
      for(int i=0; i<refsph.n; i++) dum += ssd[i];
      ssd[iter]=dum*refsph.n;
      if(iter>=2 && ssd[iter]==ssd[iter-2]) break;

      for(int s=0; s<nim; s++)
      {
         loomsk[s]=0; // leaves s out

         extract_patches(nim, refsph, im, lm_in, A, loomsk);
         detect_lm(searchsph, testsph, im[s], search_center, refsph, lm_out[s]);

         loomsk[s]=1; // add s in 
      }

      d=lm_distance(lm_in, lm_out, nim);
      lm_copy(lm_out, lm_in, nim);
      compute_landmark_cm(nim, lm_in, lmcm, loomsk);
      
      if(d==0) break;
   }

   return(d);
}

int main(int argc, char **argv)
{
   time_t time_start, time_end;

   time(&time_start);

   nifti_1_header imhdr;
   float *invT;

   char dums[DEFAULT_STRING_LENGTH];
   FILE *fp;
   char outputfile[DEFAULT_STRING_LENGTH];

   int cbd; //city block distance 

   // filename where the list of images to be registered are input 
   char imlist[DEFAULT_STRING_LENGTH]=""; 

   int nim=0; // number of images in imlist

   // An n-array of 4x4 transformation matrices which transform each of
   // the n images into the standardized PIL system
   float **TPIL;
   float **TLM;

   // An n-array of strings where the n filenames in imlist are saved
   char **imfile;
   char **lmfile;
   char **prefix;

   // temporary filename used for reading or writing specific files
   char filename[DEFAULT_STRING_LENGTH];

   int threshold=50;

   int patch_radius=DEFAULT_PATCh_RADIUS;
   int search_radius=DEFAULT_SEARCH_RADIUS;

   SHORTIM PILim[MAXIM];
   SHORTIM im[MAXIM];

   int lm_in[MAXIM][3]; // (i,j,k) coordinates of the input landmarks
   int lm_out[MAXIM][3];
   int lmcm[3]; // landmarks center of mass
   float *ssd;
   float *A; // (NxM) matrix where N=n; M=sph.n

   char loomsk[MAXIM];

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) 
      {
         case 'm':
            maxiter = atoi(optarg);
            if(maxiter<10 || maxiter>100) maxiter=MAXITER;
            break;
         case 'i':
            sprintf(imlist,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'r':
            patch_radius = atoi(optarg);
            if(patch_radius<=0 || patch_radius>MAX_RADIUS) patch_radius=DEFAULT_PATCh_RADIUS;
            break;
         case 'R':
            search_radius = atoi(optarg);
            if(search_radius<=0 || search_radius>MAX_RADIUS) search_radius=DEFAULT_SEARCH_RADIUS;
            break;
         case 't':
            threshold = atoi(optarg);
            if(threshold<0 || threshold>100) threshold=50;
            break;
         case 'o':
            sprintf(outputfile,"%s",optarg);
            opt_o=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case '?':
            print_help_and_exit();
      }
   }

   getARTHOME();

   if(argc==1) print_help_and_exit();

   if(imlist[0]=='\0')
   {
      printf("Please specify an image list using the -imlist argument\n");
      exit(0);
   }

   // find nim
   fp=fopen(imlist,"r");
   nim=0;
   while(fscanf(fp,"%s",dums) != EOF ) 
   {
         if( not_magical_nifti(dums,0)==0 ) nim++;
   }
   fclose(fp);

   if(opt_v) printf("Number of images to be registered = %d\n", nim);

   if(nim<=1) 
   {
      printf("At least two 3D T1W NIFTI images of type 'short' must be specified in %s.\n", imlist);
      exit(0);
   }

   // allocate memory for image filenames and read them into imfile
   imfile = (char **)calloc(nim, sizeof(char *));
   lmfile = (char **)calloc(nim, sizeof(char *));
   prefix = (char **)calloc(nim, sizeof(char *));
   TPIL = (float **)calloc(nim, sizeof(float *));
   TLM = (float **)calloc(nim, sizeof(float *));
   for(int i=0; i<nim; i++)
   {
      imfile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      sprintf(imfile[i],"");
      lmfile[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      prefix[i] = (char *)calloc(DEFAULT_STRING_LENGTH, sizeof(char));
      TPIL[i] = (float *)calloc(16, sizeof(float));
      TLM[i] = (float *)calloc(16, sizeof(float));
   }

   fp=fopen(imlist,"r");
   nim=0;
   while(fscanf(fp,"%s",dums) != EOF ) 
   {
         if( not_magical_nifti(dums,0)==0 ) 
         {
            strcpy(imfile[nim], dums);
            if( niftiFilename(prefix[nim], imfile[nim])==0 ) { exit(0); }
            nim++;
         }
         else
         {
            // ensure nim>0 in case the put a lm file first
            if(nim>0) strcpy(lmfile[nim-1], dums);
         }
   }
   fclose(fp);

   if(opt_v) 
   {
      printf("Input image list:\n");
      for(int i=0; i<nim; i++) 
      {
         printf("Image %d = %s\n",i+1,imfile[i]);
         if(lmfile[i][0]!='\0') printf("Image %d landmarks file = %s\n",i+1,lmfile[i]);
      }
      
   }

   /////////////////////////////////////////////////////////////////////////////////////////////
   // read PILbraincloud.nii from the $ARTHOME directory
   /////////////////////////////////////////////////////////////////////////////////////////////
   short *PILbraincloud;
   nifti_1_header PILbraincloud_hdr; 

   sprintf(filename,"%s/PILbrain.nii",ARTHOME);

   PILbraincloud = (short  *)read_nifti_image(filename, &PILbraincloud_hdr);

   if(PILbraincloud==NULL)
   {
         printf("Error reading %s, aborting ...\n", filename);
         exit(1);
   }
   /////////////////////////////////////////////////////////////////////////////////////////////
 
   if(opt_v) printf("PIL transformation:\n");
   for(int i=0; i<nim; i++)
   {
      if(opt_v) printf("Processing %s ...\n",imfile[i]);
      new_PIL_transform(imfile[i],lmfile[i],TPIL[i]);
   }

   for(int i=0; i<nim; i++)
   {
      set_dim(PILim[i], PILbraincloud_hdr);

      im[i].v = (short *)read_nifti_image(imfile[i], &imhdr);
      set_dim(im[i], imhdr);

      invT = inv4(TPIL[i]);
      PILim[i].v = resliceImage(im[i].v, im[i].nx, im[i].ny, im[i].nz, im[i].dx, im[i].dy, im[i].dz,
      PILim[i].nx, PILim[i].ny, PILim[i].nz, PILim[i].dx, PILim[i].dy, PILim[i].dz, invT, LIN);
      free(invT);

      sprintf(filename,"%s_PIL0.nii",prefix[i]);
      save_nifti_image(filename, PILim[i].v, &PILbraincloud_hdr);
   }

   SPH refsph(patch_radius);
   SPH testsph(patch_radius);
   SPH searchsph(search_radius);
   ssd = (float *)calloc(refsph.n,sizeof(float));
   A = (float *)calloc(nim*refsph.n,sizeof(float));

   //////////////////////////////////////////////////////////
   // Compute the number of seeds
   int nseeds=0;
   int *seedi, *seedj, *seedk;
   for(int i=0; i<SNX; i += delx)
   for(int j=0; j<SNY; j += dely)
   for(int k=14; k<SNZ; k += delz)
   if ( PILbraincloud[k*SNX*SNY + j*SNX + i] > threshold )
   {
      nseeds++;
   }

   seedi = (int *)calloc(nseeds, sizeof(int));
   seedj = (int *)calloc(nseeds, sizeof(int));
   seedk = (int *)calloc(nseeds, sizeof(int));

   nseeds=0;
   for(int i=0; i<SNX; i += delx)
   for(int j=0; j<SNY; j += dely)
   for(int k=14; k<SNZ; k += delz)
   if ( PILbraincloud[k*SNX*SNY + j*SNX + i] > threshold )
   {
      seedi[nseeds]=i;
      seedj[nseeds]=j;
      seedk[nseeds]=k;
      nseeds++;
   }
   //////////////////////////////////////////////////////////
  
   float **Pt;
   float **P;
   int nlm;
   int nlm_nonconvergent;
   int nlm_aligned;
   int nlm_aligned_old=0;

   Pt = (float **)calloc(nim, sizeof(float *) );
   P = (float **)calloc(nim, sizeof(float *) );

   // the Pt[m] are actually nlm*3 matrices but we allocate more rows (nseeds rows)
   // since at this time we only know that nlm<=nseeds
   for(int m=0;m<nim;m++) Pt[m] = (float *)calloc(nseeds*3, sizeof(float) );

   int sameflag;

   if(opt_v) printf("LOOC landmark detection ...\n");
   if(opt_v && maxiter!=MAXITER) printf("Maximum number of allowed iterations per seed = %d\n",maxiter);

   for(int iteration=0; iteration<maxiter; iteration++)
   {
      nlm=0;
      nlm_nonconvergent=0;
      nlm_aligned=0;

      for(int s=0; s<nseeds; s++)
      {
         for(int m=0; m<nim; m++)
         {
            lm_in[m][0]=seedi[s];
            lm_in[m][1]=seedj[s];
               lm_in[m][2]=seedk[s];
         }
   
         cbd = seek_lm(nim, refsph, PILim, lm_in, A, ssd, loomsk, lmcm, searchsph, testsph, lm_out);
   
         if(cbd>0) { nlm_nonconvergent++; continue; }
   
         /////////////////////////////////////////////////////////////////////////////////////////////////////
         // this seemed to improve things in registering 941_S_4255_B0.nii and 941_S_4255_F1.nii
         /////////////////////////////////////////////////////////////////////////////////////////////////////
         sameflag=1;
         for(int m=1; m<nim; m++)
         {
            if(lm_out[m][0] != lm_out[0][0] || lm_out[m][1] != lm_out[0][1] || lm_out[m][2] != lm_out[0][2] )
            { sameflag=0; break; }
         }
   
         //if(sameflag) { nlm_aligned++; continue;}
         if(sameflag) { nlm_aligned++;}
         /////////////////////////////////////////////////////////////////////////////////////////////////////

         //printf("\nLM %d\n",nlm);
   
         for(int m=0; m<nim; m++)
         {
   
         Pt[m][nlm*3 + 0] = lm_out[m][0]; 
         Pt[m][nlm*3 + 1] = lm_out[m][1]; 
         Pt[m][nlm*3 + 2] = lm_out[m][2]; 

            //printf("%s %d %d %d\n",imfile[m],lm_out[m][0],lm_out[m][1],lm_out[m][2]);
         }
   
         nlm++;
      }

      if(nlm_aligned <= nlm_aligned_old && iteration>0) break;
      else nlm_aligned_old=nlm_aligned;

      if(opt_v) printf("Iteration %d ...\n", iteration+1);
      //if(opt_v) printf("Number of seeds = %d\n", nseeds);
      if(opt_v) printf("Number of noncovergent seeds = %d\n",nlm_nonconvergent);
      if(opt_v) printf("Number of aligned landmarks = %d\n",nlm_aligned);
      if(opt_v) printf("Number of non-aligned landmarks = %d\n",nlm-nlm_aligned);

      if(nlm<6)
      {
         printf("Warning: insufficient LOOC detected. Skipping the generalized Procrustes anlaysis ...\n");
         break;
      }
      else
      {
         float *Q;
         float dum[3];
         float *Ptmp;
         float *Qtmp;
   
         if(opt_v) printf("Generalized Procrustes anlaysis ...\n");
   
         for(int m=0; m<nim; m++)
         {
            P[m] = (float *)calloc(3*nlm, sizeof(float) );
            transpose_matrix(Pt[m], nlm,  3, P[m]);
            convert_to_xyz(P[m], nlm, PILim[m]);
         }
   
         Q = (float *)calloc(3*nlm, sizeof(float) );
         Qtmp = (float *)calloc(3*nlm, sizeof(float) );
         Ptmp = (float *)calloc(3*nlm, sizeof(float) );
   
         // compute initial Q
         for(int k=0; k<3*nlm; k++)
         {
               Q[k]=0.0;
               for(int m=0; m<nim; m++) Q[k] += P[m][k];
               Q[k]/=nim;
         }
   
         for(int i=0; i<100; i++)
         {
            for(int m=0; m<nim; m++) 
            {
               for(int k=0; k<3*nlm; k++) { Ptmp[k]=P[m][k]; Qtmp[k]=Q[k]; }
   
               Procrustes(Qtmp, nlm, Ptmp, TLM[m]);
   
   //            printMatrix(TLM[m],4,4,"",NULL);
            }
   
            // update Q: find a new Q as the average of transformed P[m]
            for(int k=0; k<nlm; k++) 
            {
               Q[k]=Q[k+nlm]=Q[k+nlm*2]=0.0;
               for(int m=0; m<nim; m++) 
               {
                  dum[0] = P[m][k];
                  dum[1] = P[m][nlm + k];
                  dum[2] = P[m][2*nlm + k];
   
                  Q[k]         += dum[0]*TLM[m][0] + dum[1]*TLM[m][1] + dum[2]*TLM[m][2]  + TLM[m][3];
                  Q[nlm + k]   += dum[0]*TLM[m][4] + dum[1]*TLM[m][5] + dum[2]*TLM[m][6]  + TLM[m][7];
                  Q[2*nlm + k] += dum[0]*TLM[m][8] + dum[1]*TLM[m][9] + dum[2]*TLM[m][10] + TLM[m][11];
               }
               Q[k]       /= nim;
               Q[k+nlm]   /= nim;
               Q[k+2*nlm] /= nim;
            }
         }
         free(Q);
         free(Ptmp);
         free(Qtmp);
      }
   
      // update TPIL 
      for(int i=0; i<nim; i++)
      {
         multi(TLM[i],4,4,TPIL[i],4,4,TPIL[i]);
      }

      // update PILim
      for(int i=0; i<nim; i++)
      {
         free(PILim[i].v);
         invT = inv4(TPIL[i]);
         PILim[i].v = resliceImage(im[i].v, im[i].nx, im[i].ny, im[i].nz, im[i].dx, im[i].dy, im[i].dz,
         PILim[i].nx, PILim[i].ny, PILim[i].nz, PILim[i].dx, PILim[i].dy, PILim[i].dz, invT, LIN);
         free(invT);
      }
   }
   // end of iterations
   
#if 0
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   fp=fopen("secret","w");
   srand48(time(NULL));
   double random_number;
   float random_mat[16];
   random_number=drand48();
   fprintf(fp,"random_number 1 = %lf\n",random_number);
   set_to_I(random_mat,4);
   if(random_number<=(1./6.)) 
   {
      random_mat[3]=0.5;
   }
   else if(random_number<=(2./6.)) 
   {
      random_mat[7]=0.5;
   }
   else if(random_number<=(3./6.)) 
   {
      random_mat[11]=0.5;
   }
   else if(random_number<=(4./6.)) 
   {
      rotate(random_mat, .5*3.14159/180., 1.0, 0.0, 0.0);
   }
   else if(random_number<=(5./6.)) 
   {
      rotate(random_mat, .5*3.14159/180., 0.0, 1.0, 0.0);
   }
   else
   {
      rotate(random_mat, .5*3.14159/180., 0.0, 0.0, 1.0);
   }

   random_number=drand48();
   fprintf(fp,"random_number 2 = %lf\n",random_number);
   if(random_number<=0.25)
   {
      multi(random_mat,4,4,TPIL[2],4,4,TPIL[2]);
      free(PILim[2].v);
      invT = inv4(TPIL[2]);
      PILim[2].v = resliceImage(im[2].v, im[2].nx, im[2].ny, im[2].nz, im[2].dx, im[2].dy, im[2].dz,
      PILim[2].nx, PILim[2].ny, PILim[2].nz, PILim[2].dx, PILim[2].dy, PILim[2].dz, invT, LIN);
      free(invT);
      fprintf(fp,"F0\n");
   }
   else if(random_number<=0.50)
   {
      multi(random_mat,4,4,TPIL[3],4,4,TPIL[3]);
      free(PILim[3].v);
      invT = inv4(TPIL[3]);
      PILim[3].v = resliceImage(im[3].v, im[3].nx, im[3].ny, im[3].nz, im[3].dx, im[3].dy, im[3].dz,
      PILim[3].nx, PILim[3].ny, PILim[3].nz, PILim[3].dx, PILim[3].dy, PILim[3].dz, invT, LIN);
      free(invT);
      fprintf(fp,"F1\n");
   }
   else if(random_number<=0.75)
   {
      multi(random_mat,4,4,TPIL[1],4,4,TPIL[1]);
      free(PILim[1].v);
      invT = inv4(TPIL[1]);
      PILim[1].v = resliceImage(im[1].v, im[1].nx, im[1].ny, im[1].nz, im[1].dx, im[1].dy, im[1].dz,
      PILim[1].nx, PILim[1].ny, PILim[1].nz, PILim[1].dx, PILim[1].dy, PILim[1].dz, invT, LIN);
      free(invT);
      fprintf(fp,"B1\n");
   }
   else
   {
      multi(random_mat,4,4,TPIL[0],4,4,TPIL[0]);
      free(PILim[0].v);
      invT = inv4(TPIL[0]);
      PILim[0].v = resliceImage(im[0].v, im[0].nx, im[0].ny, im[0].nz, im[0].dx, im[0].dy, im[0].dz,
      PILim[0].nx, PILim[0].ny, PILim[0].nz, PILim[0].dx, PILim[0].dy, PILim[0].dz, invT, LIN);
      free(invT);
      fprintf(fp,"B0\n");
   }
   fclose(fp);
   ///////////////////////////////////////////////////////////////////////////////////////////////////
#endif

   // save PIL images
   for(int i=0; i<nim; i++)
   {
      sprintf(filename,"%s_PIL.nii",prefix[i]);
      save_nifti_image(filename, PILim[i].v, &PILbraincloud_hdr);
   }

   // save transformation matrix
   for(int i=0; i<nim; i++)
   {
      sprintf(filename,"%s_PIL.mrx",prefix[i]);
      fp=fopen(filename,"w");
      printMatrix(TPIL[i],4,4,"",fp);
      fclose(fp);
   }

   // free memory
   for(int i=0; i<nim; i++)
   {
      free(imfile[i]);
      free(prefix[i]);
      free(TPIL[i]);
      free(TLM[i]);
      free(PILim[i].v);
   }
   free(ssd);
   free(PILbraincloud);

   time(&time_end);
   if(opt_v) printf("Execution time = %d minutes and %d seconds\n",(time_end-time_start)/60,
   (time_end-time_start)%60);

   return 0;
}
