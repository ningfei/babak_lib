// The mail update in this version of the program is that it utilizes multiple processors using MPI
// Added the -cc option

#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <unistd.h>
#include <volume.h>
#include <spm_analyze.h>
#include <nifti1_io.h>
#include <babak_lib.h>
#include <stats.h>
#include <minmax.h>
#include <smooth.h>
#include <mpi.h>

#define YES 1
#define NO 0

#define NX 512
#define NY 512
#define NP 512*512

#define UPPER_LEFT_i 140
#define UPPER_LEFT_j 160
#define LOWER_RIGHT_i 360
#define LOWER_RIGHT_j 280

//////////////////////////////////////////////////////////////////////////////////////////////////
float AC[4]={0.0, 0.0, 0.0, 1.0};
float PC[4]={0.0, 0.0, 0.0, 1.0};
float VSPS[4]={0.0, 0.0, 0.0, 1.0};

//int np;
float dx, dy, dz;
int N2;	// N*N

float *Xwarp, *Ywarp, *Zwarp;
int niter=4;

float *ARobj;	// array extracted from object and target images
float *ARtrg;

int Wx, Wy;
int Lx, Ly;

nifti_1_header sub_hdr;
nifti_1_header output_hdr; // root only
short *subject_volume;
int Snx, Sny, Snz;
float Sdx, Sdy, Sdz;

float Tacpc[16];

float hampel_origin[2];
float hampel_axis[2];
int inferior_genu[2];
int inferior_splenium[2];
int anterior_point[2];
int posterior_point[2];
int posterior_genu[2];
int rostrum[2];
int cc_length_fifth;
int cc_length_third;
int cc_length_half;
int cc_length;

int *cci, *ccj;

int *medi, *medj; // (i,j) components of the medial axis pixels 0, 1, ..., nm-1
int nm=0;  // number of pixels on the medial axis
short *dist;

// normal unit vectors to the medial axis
// program is designed such that the normal vectors to the medial axis always point towards 
// the upper boundary of the CC.
float *normi, *normj; 

int *ubi, *ubj, *lbi, *lbj;
int nb; // number of border pixels
int nub; // number of upper boundary pixels
int nlb; // number of lower boundary pixels

float *lbindx, *ubindx;

int acpoint;
int pcpoint;

float ACx;
float PCx;

float ACi;
float PCi;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-A", 1,  'L'},
   {"-version", 0,  'V'},
   {"-h", 0,  'h'},
   {"-help", 0,  'h'},
   {"-v", 0,  'v'},
   {"-verbose", 0,  'v'},
   {"-a", 1, 'a'},
   {"-atlas", 1, 'a'},
   {"-o", 1, 'o'},
   {"-i", 1,  'i'},
   {"-n", 1,  'n'},
   {"-T", 1,  'T'},
   {"-csv", 1,  'c'},
   {"-mrx", 0,  'm'},
   {"-ppm", 0,  'p'},
   {"-acpc", 0,  'p'},
   {"-box", 0,  'b'},
   {"-secret", 0,  's'},
   {"-W", 0,  'W'},
   {"-Hampel", 0,  'H'},
   {"-H", 0,  'H'},
   {"-border", 0, 'B'},
   {"-cc", 1, 'C'},
   {"-VSPS",1,'S'},
   {"-AC",1,'A'},
   {"-PC",1,'P'},
   {"-t",1,'t'},
   {"-thresh",1,'t'},
   {"-threshold",1,'t'},
   {0, 0,  0}
};

int opt_a=NO; // if YES, a user-specific atlas has been given at the command line
int opt_cc=NO;
int opt_mrx=NO;
int opt_box=NO;
int opt_W=NO;
int opt_H=NO;
int opt_border=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help();

void computeWarpField(short *obj, short *trg, float sd, int bbnp, int bbnx, int bbny)
{
   int xopt, yopt;

   float Vt, Vo;
   float Sx, Sx2, Sy, Sy2, Sxy;
   float num, den;
   float CC, CCMAX;
   float *Xw, *Yw;
   float *Xtmp, *Ytmp;

   short *HRtrg;	// high res. target and object images
   short *HRobj;

   float HRdx, HRdy;
   int HRnx, HRny;

   for(int n=0; n<bbnp; n++) Xwarp[n]=Ywarp[n]=0.0;

   for(int iter=0; iter<niter; iter++)
   {
      // printf("\n\tIteration %d\n",iter);

      HRdx = (float)(dx*pow(2.0,(niter-iter-1.0)));
      HRdy = (float)(dy*pow(2.0,(niter-iter-1.0)));

      HRnx = (int)(bbnx/pow(2.0,(niter-iter-1.0)) + 0.5 );
      HRny = (int)(bbny/pow(2.0,(niter-iter-1.0)) + 0.5 );

      //printf("\tMatrix size = %d x %d (voxels)\n", HRnx, HRny);
      //printf("\tVoxel size = %8.6f x %8.6f (mm3)\n", HRdx,HRdy);

      Xw=resizeXY(Xwarp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy);
      Yw=resizeXY(Ywarp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy);

      {
         float *tmp;
         float StandDev;

         StandDev = (0.5/log(2.0)) * ( HRdx*HRdx - dx*dx );

         if(StandDev>0.0)
         {
            if( ( HRdx*HRdx - dx*dx ) > 0.0 )
            {
               StandDev=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - dx*dx ) )/dx );
            }

            tmp = smoothXY(obj,bbnx,bbny,StandDev);
            HRobj=computeReslicedImage(tmp, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
            free(tmp);
         }
         else
         {
            HRobj=computeReslicedImage(obj, bbnx, bbny, dx, dy, HRnx, HRny, HRdx, HRdy, Xw, Yw);
         }
      }

      HRtrg=resizeXY(trg, bbnx ,bbny, dx, dy, HRnx, HRny, HRdx, HRdy);

      for(int j=0; j<HRny; j++)
      for(int i=0; i<HRnx; i++)
      {
			extractArray(HRtrg, HRnx, HRny, i, j, Lx, Ly, ARtrg);

			Sy=Sy2=0.0;
			for(int n=0; n<N2; n++)
			{
				Vt = ARtrg[n];
				Sy += Vt;
				Sy2 += (Vt*Vt);
			}

			if( Sy == 0.0 )
			{
				Xw[j*HRnx + i] = 0.0;
				Yw[j*HRnx + i] = 0.0;
				continue;
			}

			CCMAX=0.0; 	// IMPORTANT: we are not interested in -tive correlations
					// if CMAX is set to -1, program given unexpected results

			xopt=yopt=0;

			for(int x=-Wx; x<=Wx; x++)
			for(int y=-Wy; y<=Wy; y++)
			{
				extractArray(HRobj, HRnx, HRny, i+x, j+y, Lx, Ly, ARobj);
	
				Sx=Sx2=Sxy=0.0;
				for(int n=0; n<N2; n++)
				{
					Vo = ARobj[n];
					Sx += Vo;
					Sx2 += (Vo*Vo);
					Sxy += (Vo*ARtrg[n]);
				}

				num = Sxy-Sx*Sy/N2;
				den = (float)sqrt( (double)(Sx2-Sx*Sx/N2)*(Sy2-Sy*Sy/N2) );
	
				if(den==0.0) continue;
	
				CC = num/den;
		
				if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; }
			}
	
			Xw[j*HRnx + i] = xopt*HRdx;
			Yw[j*HRnx + i] = yopt*HRdy;
      }
	
      Xtmp=smoothXY(Xw, HRnx, HRny, sd);
      Ytmp=smoothXY(Yw, HRnx, HRny, sd);
	
      free(Xw); free(Yw);
	
      Xw=resizeXY(Xtmp, HRnx, HRny, HRdx, HRdy, bbnx, bbny, dx, dy);
      Yw=resizeXY(Ytmp, HRnx, HRny, HRdx, HRdy, bbnx, bbny, dx, dy);
		
      free(Xtmp); free(Ytmp);
		
      for(int n=0; n<bbnp; n++)
      {
         Xwarp[n] += Xw[n];
         Ywarp[n] += Yw[n];
      }
		
      free(Xw); free(Yw);
	
      free(HRobj);
      free(HRtrg);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void update_qsform( const char *imagefilename , float *matrix)
{
   FILE *fp;
   nifti_1_header hdr; // 348 bytes
   nifti1_extender ext; // 4 bytes
   char *extension=NULL;
   int extension_size=0;
   char *data=NULL;
   int data_size=0;
   char swapflg=0;
   mat44 R;

   fp = fopen(imagefilename,"r");
   if(fp==NULL) file_open_error(imagefilename);
   fread(&hdr, sizeof(nifti_1_header), 1, fp);

   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
      swapflg=1;
   }

   //if( opt_sform)
   {
      hdr.sform_code = NIFTI_XFORM_TALAIRACH;
      hdr.srow_x[0]=matrix[0]; hdr.srow_x[1]=matrix[1]; hdr.srow_x[2]=matrix[2]; hdr.srow_x[3]=matrix[3];
      hdr.srow_y[0]=matrix[4]; hdr.srow_y[1]=matrix[5]; hdr.srow_y[2]=matrix[6]; hdr.srow_y[3]=matrix[7];
      hdr.srow_z[0]=matrix[8]; hdr.srow_z[1]=matrix[9]; hdr.srow_z[2]=matrix[10]; hdr.srow_z[3]=matrix[11];
   }

   //if( opt_qform)
   {
      hdr.qform_code = NIFTI_XFORM_TALAIRACH;
      R.m[0][0]=matrix[0];  R.m[0][1]=matrix[1];  R.m[0][2]=matrix[2];  R.m[0][3]=matrix[3];
      R.m[1][0]=matrix[4];  R.m[1][1]=matrix[5];  R.m[1][2]=matrix[6];  R.m[1][3]=matrix[7];
      R.m[2][0]=matrix[8];  R.m[2][1]=matrix[9];  R.m[2][2]=matrix[10]; R.m[2][3]=matrix[11];
      R.m[3][0]=matrix[12]; R.m[3][1]=matrix[13]; R.m[3][2]=matrix[14]; R.m[3][3]=matrix[15];

      nifti_mat44_to_quatern( R,  &(hdr.quatern_b), &(hdr.quatern_c), &(hdr.quatern_d),
      &(hdr.qoffset_x), &(hdr.qoffset_y), &(hdr.qoffset_z), 
      &(hdr.pixdim[1]), &(hdr.pixdim[2]), &(hdr.pixdim[3]), &(hdr.pixdim[0]));
   }

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      fread(&ext, sizeof(nifti1_extender), 1, fp);

      extension_size = (int)(hdr.vox_offset)-352;

      if( extension_size > 0 )
      {
         extension = (char *)calloc(extension_size, 1);
         fread(extension, 1, extension_size, fp);
      }

      data_size = 1;
      for(int i=1; i<=hdr.dim[0]; i++)
      {
         data_size *= hdr.dim[i];
      }
      data_size *= (hdr.bitpix/8);

      if( data_size > 0 )
      {
         data = (char *)calloc(data_size, 1);
         fread(data, 1, data_size, fp);
      }
   }

   fclose(fp);

   fp = fopen(imagefilename,"w");
   if(fp==NULL) file_open_error(imagefilename);

   if(swapflg)
   {
      swapniftiheader(&hdr);
   }

   fwrite(&hdr, sizeof(nifti_1_header), 1, fp);

   if( hdr.magic[0]=='n' && hdr.magic[1]=='+' && hdr.magic[2]=='1' )
   {
      fwrite(&ext, sizeof(nifti1_extender), 1, fp);
      if( extension_size > 0 )
      {
         fwrite(extension, 1, extension_size, fp);
         delete extension;
      }
      if( data_size > 0 )
      {
         fwrite(data, 1, data_size, fp);
         delete data;
      }
   }
   fclose(fp);
}
//////////////////////////////////////////////////////////////////////////////////////////////////

// this function should only be run under procID=root
short *find_subject_msp(char *imagefilename, char *prefix)
{
   DIM input_dim, output_dim;
   short *msp;
   char orientation[4]="";
   char modelfile[1024]="";
   int opt_T2=NO;

   float Tmsp[16]; // transforms volOrig to MSP aligned PIL orientation

   float ac[4], pc[4];  

   float *invT;

   char outputfilename[512];
   FILE *fp;

   detect_AC_PC_MSP(imagefilename, orientation, modelfile, AC, PC, VSPS, Tmsp, 0, opt_T2);

   input_dim.nx = Snx;
   input_dim.ny = Sny;
   input_dim.nz = Snz;
   input_dim.dx = Sdx;
   input_dim.dy = Sdy;
   input_dim.dz = Sdz;

   output_dim.nx = NX;
   output_dim.ny = NY;
   output_dim.nz = 1;
   output_dim.nt = 1;
   output_dim.dx = dx;
   output_dim.dy = dy;
   output_dim.dz = dz;
   output_dim.dt = 0.0;

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   for(int i=0; i<4; i++) ac[i] = AC[i];
   for(int i=0; i<4; i++) pc[i] = PC[i];

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   orig_ijk_to_pil_xyz(Tmsp, input_dim, ac, pc);
   ACPCtransform(Tacpc, Tmsp, ac, pc, 0);

   // this part of the code located the ac and pc locations on the MSP-AC-PC aligned sagittal slice in (i,j) coordinates.
   // the actuall ac location is actuall approximately acpoint+0.5.  also pc is approx. pcpoint+0.5 
   // on the output PPM image, we mark points acpoint, acpoint+1, pcpoint and pcpoint+1 on rows 255 and 256.
   {
      float I2X[16];
      float X2I[16];

      for(int i=0; i<4; i++) ac[i] = AC[i];
      for(int i=0; i<4; i++) pc[i] = PC[i];

      ijk2xyz(I2X, input_dim.nx, input_dim.ny, input_dim.nz, input_dim.dx, input_dim.dy, input_dim.dz);

      multi(I2X, 4, 4,  ac, 4,  1, ac);
      multi(I2X, 4, 4,  pc, 4,  1, pc);

      multi(Tacpc, 4, 4,  ac, 4,  1, ac);
      multi(Tacpc, 4, 4,  pc, 4,  1, pc);

      ACx=ac[0]; 
      PCx=pc[0];

      xyz2ijk(X2I, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);

      multi(X2I, 4, 4,  ac, 4,  1, ac);
      multi(X2I, 4, 4,  pc, 4,  1, pc);

      acpoint = (int)ac[0];
      pcpoint = (int)pc[0];

      ACi=ac[0];
      PCi=pc[0];
   }

   invT = inv4(Tacpc);
   msp = resliceImage(subject_volume,input_dim, output_dim,invT,LIN);
   free(invT);
   
   // If -mrx option is selected, save the matrix that makes the input image to MSP/AC-PC aligned
   if(opt_mrx)
   {
      sprintf(outputfilename,"%s_msp.mrx",prefix);
      fp = fopen(outputfilename,"w");

      if(fp != NULL)
      {
         printMatrix(Tacpc, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
         fclose(fp);
      }
      else
      {
         printf("Cound not write to %s.\n", outputfilename);
      }
   }

   {
      float T_ijk2xyz[16];
      float PIL2RAS[16];

      output_hdr = read_NIFTI_hdr(imagefilename);
      output_hdr.pixdim[1]=output_dim.dx; 
      output_hdr.pixdim[2]=output_dim.dy; 
      output_hdr.pixdim[3]=output_dim.dz;
      output_hdr.dim[1]=output_dim.nx; 
      output_hdr.dim[2]=output_dim.ny; 
      output_hdr.dim[3]=output_dim.nz;
      output_hdr.magic[0]='n'; output_hdr.magic[1]='+'; output_hdr.magic[2]='1';
      sprintf(output_hdr.descrip,"Created by ART yuki");

      sprintf(outputfilename,"%s_msp.nii",prefix);
      save_nifti_image(outputfilename, msp, &output_hdr);

      //////////////////////////////////////////////////////////////////////////////////
      // This part of the code adjusts the SFORM matrix of the output image

      inversePILtransform("RAS", PIL2RAS);

      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      update_qsform( (const char *)outputfilename, Tacpc );
      //////////////////////////////////////////////////////////////////////////////////
   }

   return(msp);
}

// this function should only be run under procID=root
short *find_subject_msp(char *imagefilename, char *prefix, char *msp_transformation_file)
{
   DIM input_dim, output_dim;
   short *msp;
   char orientation[4]="";
   char modelfile[1024]="";
   int opt_T2=NO;

   float Tmsp[16]; // transforms volOrig to MSP aligned PIL orientation

   float ac[4], pc[4];  

   float *invT;

   char outputfilename[512];

   detect_AC_PC_MSP(imagefilename, orientation, modelfile, AC, PC, VSPS, Tmsp, 0, opt_T2);

   input_dim.nx = Snx;
   input_dim.ny = Sny;
   input_dim.nz = Snz;
   input_dim.dx = Sdx;
   input_dim.dy = Sdy;
   input_dim.dz = Sdz;

   output_dim.nx = NX;
   output_dim.ny = NY;
   output_dim.nz = 1;
   output_dim.nt = 1;
   output_dim.dx = dx;
   output_dim.dy = dy;
   output_dim.dz = dz;
   output_dim.dt = 0.0;

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   for(int i=0; i<4; i++) ac[i] = AC[i];
   for(int i=0; i<4; i++) pc[i] = PC[i];

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   orig_ijk_to_pil_xyz(Tmsp, input_dim, ac, pc);
   ACPCtransform(Tacpc, Tmsp, ac, pc, 0);

   // this part of the code locates the ac and pc locations on the MSP-AC-PC aligned sagittal slice in (i,j) coordinates.
   // the actuall ac location is actually approximately acpoint+0.5.  also pc is approx. pcpoint+0.5 
   // on the output PPM image, we mark points acpoint, acpoint+1, pcpoint and pcpoint+1 on rows 255 and 256.
   {
      float I2X[16];
      float X2I[16];

      for(int i=0; i<4; i++) ac[i] = AC[i];
      for(int i=0; i<4; i++) pc[i] = PC[i];

      ijk2xyz(I2X, input_dim.nx, input_dim.ny, input_dim.nz, input_dim.dx, input_dim.dy, input_dim.dz);

      multi(I2X, 4, 4,  ac, 4,  1, ac);
      multi(I2X, 4, 4,  pc, 4,  1, pc);

      multi(Tacpc, 4, 4,  ac, 4,  1, ac);
      multi(Tacpc, 4, 4,  pc, 4,  1, pc);

      ACx=ac[0]; 
      PCx=pc[0];

      xyz2ijk(X2I, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);

      multi(X2I, 4, 4,  ac, 4,  1, ac);
      multi(X2I, 4, 4,  pc, 4,  1, pc);

      acpoint = (int)ac[0];
      pcpoint = (int)pc[0];

      ACi=ac[0];
      PCi=pc[0];
   }

   loadTransformation(msp_transformation_file, Tacpc);

   invT = inv4(Tacpc);
   msp = resliceImage(subject_volume,input_dim, output_dim,invT,LIN);
   free(invT);

   {
      float T_ijk2xyz[16];
      float PIL2RAS[16];

      output_hdr = read_NIFTI_hdr(imagefilename);
      output_hdr.pixdim[1]=output_dim.dx; 
      output_hdr.pixdim[2]=output_dim.dy; 
      output_hdr.pixdim[3]=output_dim.dz;
      output_hdr.dim[1]=output_dim.nx; 
      output_hdr.dim[2]=output_dim.ny; 
      output_hdr.dim[3]=output_dim.nz;
      output_hdr.magic[0]='n'; output_hdr.magic[1]='+'; output_hdr.magic[2]='1';
      sprintf(output_hdr.descrip,"Created by ART yuki");

      sprintf(outputfilename,"%s_msp.nii",prefix);
      save_nifti_image(outputfilename, msp, &output_hdr);

      //////////////////////////////////////////////////////////////////////////////////
      // This part of the code adjusts the SFORM matrix of the output image

      inversePILtransform("RAS", PIL2RAS);

      ijk2xyz(T_ijk2xyz, output_dim.nx, output_dim.ny, output_dim.nz, output_dim.dx, output_dim.dy, output_dim.dz);
      multi(PIL2RAS, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      update_qsform( (const char *)outputfilename, Tacpc );
      //////////////////////////////////////////////////////////////////////////////////
   }

   return(msp);
}

//////////////////////////////////////////////////////////////////////////////////

void output_ppm(short *trg, short *cc_est, const char *prefix) 
{
   char outputfile[1024]="";
   short min=0, max=0;

   ////////////////////////////////////////////////////////////
   unsigned char *R, *G, *B;

   R = (unsigned char *)calloc(NP, 1);
   if(R == NULL)
   {
      memory_allocation_error("R");
   }

   G = (unsigned char *)calloc(NP, 1);
   if(G == NULL)
   {
      memory_allocation_error("G");
   }

   B = (unsigned char *)calloc(NP, 1);
   if(B == NULL)
   {
      memory_allocation_error("B");
   }

   minmax(trg, NP, min, max);

   //////////////////////////////////////////////////////////////////////////////////
   // outputs visualization of Witelson's subdivisions
   if(opt_W)
   {
      for(int i=0; i<NX; i++)
      for(int j=0; j<NY; j++)
      {
         int v;

         v = j*NX + i;

         if(cc_est[v]>0)
         {
            if(i >= (posterior_point[0]-cc_length_fifth) ) // W7 
            {
               R[v] = 143;
               G[v] = 0;
               B[v] = 255;
            }
            else if(i >= (posterior_point[0]-cc_length_third) ) // W6 
            {
               R[v] = 75;
               G[v] = 0;
               B[v] = 130;
            }
            else if(i >= (posterior_point[0]-cc_length_half) ) // W5 
            {
               R[v] = 0;
               G[v] = 0;
               B[v] = 255;
            }
            else if(i >= (anterior_point[0]+cc_length_third) ) // W4 
            {
               R[v] = 0;
               G[v] = 255;
               B[v] = 0;
            } 
            else if(i <= posterior_genu[0] ) // W2 
            {
               R[v] = 255;
               G[v] = 127;
               B[v] = 0;
            } 
            else if(j >= posterior_genu[1] ) // W1 
            {
               R[v] = 255;
               G[v] = 0;
               B[v] = 0;
            } 
            else // w3 
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 0;
            }
         }
         else
         {
            R[v] = (unsigned char)(trg[v]*255.0/max);
            G[v] = (unsigned char)(trg[v]*255.0/max);
            B[v] = (unsigned char)(trg[v]*255.0/max);
         }
      }

      sprintf(outputfile,"%s_cc_witelson.ppm",prefix);
      save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }
   //////////////////////////////////////////////////////////////////////////////////
   
   //////////////////////////////////////////////////////////////////////////////////
   // outputs visualization of Hampel's subdivisions
   if(opt_H)
   {
      float d;
      float pi;
      float i0, j0;
      float theta, costheta;
      pi = 4.0*atanf(1.0);

      i0 = hampel_origin[0];
      j0 = hampel_origin[1];

      for(int i=0; i<NX; i++)
      for(int j=0; j<NY; j++)
      {
         int v;

         v = j*NX + i;

         if(cc_est[v]>0)
         {

            d = sqrtf ( (i-i0)*(i-i0) + (j-j0)*(j-j0) );
            costheta = (i-i0)*hampel_axis[0]/d + (j-j0)*hampel_axis[1]/d ;
            if(costheta > 1.0 ) costheta=1.0;
            if(costheta < -1.0 ) costheta=-1.0;
            theta = acosf( costheta );

            if( theta <= pi/5.0) // C5 (color from line 37 of rainbow)
            {
               R[v] = 0;
               G[v] = 0;
               B[v] = 255;
            }
            else if(theta <= 2.0*pi/5.0 ) // C4 (color from line 74 of rainbow)
            {
               R[v] = 0;
               G[v] = 255;
               B[v] = 0;
            }
            else if(theta <= 3.0*pi/5.0 ) // C3 (color from line 74 of rainbow)
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 0;
            }
            else if(theta <= 4.0*pi/5.0 ) // C2 (color from line 74 of rainbow)
            {
               R[v] = 255;
               G[v] = 127;
               B[v] = 0;
            } 
            else if(theta <= pi)  // C1 (color from line 185 of rainbow)
            {
               R[v] = 255;
               G[v] = 0;
               B[v] = 0;
            }
            else
            {
               R[v] = 255;
               G[v] = 255;
               B[v] = 255;
            }
         }
         else
         {
            R[v] = (unsigned char)(trg[v]*255.0/max);
            G[v] = (unsigned char)(trg[v]*255.0/max);
            B[v] = (unsigned char)(trg[v]*255.0/max);
         }
      }

      int O[2];
      O[0] = (int)rintf( hampel_origin[0] );
      O[1] = (int)rintf( hampel_origin[1] );
      for(int i=O[0]-3; i<=O[0]+3; i++)
      {
         R[O[1]*NX+i] = 0;
         G[O[1]*NX+i] = 0;
         B[O[1]*NX+i] = 255;
      }

      for(int j=O[1]-3; j<=O[1]+3; j++)
      {
         R[j*NX+O[0]] = 0;
         G[j*NX+O[0]] = 0;
         B[j*NX+O[0]] = 255;
      }

      sprintf(outputfile,"%s_cc_hampel.ppm",prefix);
      save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }
   //////////////////////////////////////////////////////////////////////////////////

   for(int v=0; v<NP; v++)
   {
      R[v] = (unsigned char)(trg[v]*255.0/max);
      G[v] = (unsigned char)(trg[v]*255.0/max);
      B[v] = (unsigned char)(trg[v]*255.0/max);
   }

   for(int i=1; i<NX-1; i++)
   for(int j=1; j<NY-1; j++)
   {
      if(cc_est[j*NX + i]>0 && (cc_est[j*NX+i-1]==0 || cc_est[j*NX+i+1]==0 || cc_est[(j-1)*NX+i]==0 || cc_est[(j+1)*NX+i]==0)  )
      {
         R[j*NX+i] = 255;
         G[j*NX+i] = 0;
         B[j*NX+i] = 0;
      }
   }

//
// The following two lines if enables saves the MSP image with a red CC border
//
   if(opt_border)
   {
      sprintf(outputfile,"%s_cc_border.ppm",prefix);
      save_as_ppm((const char *)outputfile, NX, NY, R, G, B);
   }

   free(R);
   free(G);
   free(B);

   return;
}

//////////////////////////////////////////////////////////////////////////////////

void output_bounding_box_ppm(short *trg, short *cc, const char *prefix) 
{
   char outputfile[1024]="";
//   short min, max;
   int min=0, max=0;

   ////////////////////////////////////////////////////////////
   unsigned char *R, *G, *B;

   R = (unsigned char *)calloc(NP, 1);
   if(R == NULL)
   {
      memory_allocation_error("R");
   }

   G = (unsigned char *)calloc(NP, 1);
   if(G == NULL)
   {
      memory_allocation_error("G");
   }

   B = (unsigned char *)calloc(NP, 1);
   if(B == NULL)
   {
      memory_allocation_error("B");
   }

   // minmax(trg, NP, &min, &max);
   setLowHigh(trg, NP, &min, &max);

   for(int v=0; v<NP; v++)
   {
      if(trg[v]>max) 
      {
         R[v] = (unsigned char)(255.0);
         G[v] = (unsigned char)(255.0);
         B[v] = (unsigned char)(255.0);
      }
      else
      {
         R[v] = (unsigned char)(trg[v]*255.0/max);
         G[v] = (unsigned char)(trg[v]*255.0/max);
         B[v] = (unsigned char)(trg[v]*255.0/max);
      }
   }


   //////////////////////////////////////////////////////////////////////////////
   // draws the boudning box in *cc.ppm image
   //////////////////////////////////////////////////////////////////////////////
   if(opt_box)
   {
      for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
      {
         for(int j=UPPER_LEFT_j; j<=UPPER_LEFT_j+1; j++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }

         for(int j=LOWER_RIGHT_j-1; j<=LOWER_RIGHT_j; j++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }
      }

      for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
      {
         for(int i=UPPER_LEFT_i; i<=UPPER_LEFT_i+1; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }

         for(int i=LOWER_RIGHT_i-1; i<=LOWER_RIGHT_i; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 255;
            B[j*NX+i] = 255;
         }
      }
   }
   /////////////////////////////////////////////////////////////////////

   if( acpoint<5 || acpoint>(NX-5)) acpoint = NX/2;
   if( pcpoint<5 || pcpoint>(NX-5)) pcpoint = NX/2;
   // mark the AC and PC locations
   for(int j=(NY/2-1)-5; j<=(NY/2)+5; j++)
   {

      if(j==(NY/2-1) || j==(NY/2))
      {
         for(int i=acpoint-5; i<=acpoint+1+5; i++)
         {
            R[j*NX+i] = 0;
            G[j*NX+i] = 255;
            B[j*NX+i] = 0;
         }

         for(int i=pcpoint-5; i<=pcpoint+1+5; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 0;
            B[j*NX+i] = 0;
         }
      }
      else
      {
         for(int i=acpoint; i<=acpoint+1; i++)
         {
            R[j*NX+i] = 0;
            G[j*NX+i] = 255;
            B[j*NX+i] = 0;
         }

         for(int i=pcpoint; i<=pcpoint+1; i++)
         {
            R[j*NX+i] = 255;
            G[j*NX+i] = 0;
            B[j*NX+i] = 0;
         }
      }
   }

 
   /////////////////////////////////////////////////////////////////////////////////////
   // mark the CC border 
   // upper and lower boundaries are colored differently
   /////////////////////////////////////////////////////////////////////////////////////
   for(int b=0; b<nub; b++)
   {
      R[ ubj[b]*NX+ubi[b]] = 255;
      G[ ubj[b]*NX+ubi[b]] = 255;
      B[ ubj[b]*NX+ubi[b]] = 0;
   }
   for(int b=0; b<nlb; b++)
   {
      R[ lbj[b]*NX+lbi[b]] = 0;
      G[ lbj[b]*NX+lbi[b]] = 255;
      B[ lbj[b]*NX+lbi[b]] = 255;
   }

   /////////////////////////////////////////////////////////////////////////////////////
   // Mark the medial axis as red

   for(int v=0; v<nm; v++)
   {
      R[ medj[v]*NX + medi[v]] = 255;
      G[ medj[v]*NX + medi[v]] = 0;
      B[ medj[v]*NX + medi[v]] = 0;
   }

   //////////////////////////////////////////////////////////////////////
   // mark rostrum, most anterior and most posterior points etc.
   
   for(int i=rostrum[0]-3; i<=rostrum[0]+3; i++)
   {
      R[rostrum[1]*NX+i] = 0;
      G[rostrum[1]*NX+i] = 0;
      B[rostrum[1]*NX+i] = 255;
   }

   for(int j=rostrum[1]-3; j<=rostrum[1]+3; j++)
   {
      R[j*NX+rostrum[0]] = 0;
      G[j*NX+rostrum[0]] = 0;
      B[j*NX+rostrum[0]] = 255;
   }
   
   for(int i=anterior_point[0]-3; i<=anterior_point[0]+3; i++)
   {
      R[anterior_point[1]*NX+i] = 0;
      G[anterior_point[1]*NX+i] = 0;
      B[anterior_point[1]*NX+i] = 255;
   }

   for(int j=anterior_point[1]-3; j<=anterior_point[1]+3; j++)
   {
      R[j*NX+anterior_point[0]] = 0;
      G[j*NX+anterior_point[0]] = 0;
      B[j*NX+anterior_point[0]] = 255;
   }

   for(int i=posterior_point[0]-3; i<=posterior_point[0]+3; i++)
   {
      R[posterior_point[1]*NX+i] = 0;
      G[posterior_point[1]*NX+i] = 0;
      B[posterior_point[1]*NX+i] = 255;
   }

   for(int j=posterior_point[1]-3; j<=posterior_point[1]+3; j++)
   {
      R[j*NX+posterior_point[0]] = 0;
      G[j*NX+posterior_point[0]] = 0;
      B[j*NX+posterior_point[0]] = 255;
   }

   for(int i=inferior_genu[0]-3; i<=inferior_genu[0]+3; i++)
   {
      R[inferior_genu[1]*NX+i] = 0;
      G[inferior_genu[1]*NX+i] = 0;
      B[inferior_genu[1]*NX+i] = 255;
   }

   for(int j=inferior_genu[1]-3; j<=inferior_genu[1]+3; j++)
   {
      R[j*NX+inferior_genu[0]] = 0;
      G[j*NX+inferior_genu[0]] = 0;
      B[j*NX+inferior_genu[0]] = 255;
   }

   for(int i=inferior_splenium[0]-3; i<=inferior_splenium[0]+3; i++)
   {
      R[inferior_splenium[1]*NX+i] = 0;
      G[inferior_splenium[1]*NX+i] = 0;
      B[inferior_splenium[1]*NX+i] = 255;
   }

   for(int j=inferior_splenium[1]-3; j<=inferior_splenium[1]+3; j++)
   {
      R[j*NX+inferior_splenium[0]] = 0;
      G[j*NX+inferior_splenium[0]] = 0;
      B[j*NX+inferior_splenium[0]] = 255;
   }

   for(int i=posterior_genu[0]-3; i<=posterior_genu[0]+3; i++)
   {
      R[posterior_genu[1]*NX+i] = 0;
      G[posterior_genu[1]*NX+i] = 0;
      B[posterior_genu[1]*NX+i] = 255;
   }

   for(int j=posterior_genu[1]-3; j<=posterior_genu[1]+3; j++)
   {
      R[j*NX+posterior_genu[0]] = 0;
      G[j*NX+posterior_genu[0]] = 0;
      B[j*NX+posterior_genu[0]] = 255;
   }
   //////////////////////////////////////////////////////////////////////

   sprintf(outputfile,"%s_cc.ppm",prefix);
   save_as_ppm((const char *)outputfile, NX, NY, R, G, B);

   free(R);
   free(G);
   free(B);

   return;
}

//////////////////////////////////////////////////////////////////////////////////

// NOTE: The commented out lines reflect secret options
void print_help()
{
   printf("\nUsage: yuki [-v -version -h -o <output prefix> -csv <csvfile> -Hampel -W -acpc -cc] -i <subject volume>\n\n"

   "-v : enables verbose mode\n\n"

   "-version : reports version\n\n"

   "-h : prints help message\n\n"

   "-o <output prefix> : prefix for naming output files (default: same as input prefix)\n\n"

   "-csv <csvfile> : CC measurments (area, perimeter, etc.) will be appended to this file in\n"
   "comma-separated values (CSV) format (default: <output prefix>.csv)\n\n"

   "-Hampel: outputs <prefix>_cc_hampel.ppm as well as 5 Hampel sub-areas\n\n"

   "-W: outputs <prefix>_cc_witelson.ppm as well as 7 Witelson sub-areas\n\n"

   "-acpc: outputs files <prefix>_ACPC.txt, <prefix>_ACPC_axial.ppm, and <prefix>_ACPC_sagittal.ppm\n"
   "that show the results of AC/PC and MSP detection\n\n"

   "-cc <corrected_cc.nii>: This option is used when the out binary CC image is corrected manually and we\n"
   "need to recalculate the CC related measurements (area, circularlity, etc.) for the corrected image\n\n"

   "-i <mprage>.nii : the 3D MRI volume (short int NIFTI format) on which the corpus callosum\n"
   "is to be located\n\n");

   return;
}

void print_secret_help()
{
   printf("\nUsage: yuki [-version -h -o <output prefix> -csv <csvfile>] -i <subject volume>\n\n"
   "[-v -mrx -ppm -box -W -border -a <atlas> -n <# atlases>]\n"
   "-version : reports version\n\n"
   "-v : enables verbose mode\n\n"
   "-mrx : saves the transformation matrix that makes the input image MSP/AC-PC aligned\n\n"
   "-ppm (same as -acpc): saves files related to AC/PC and MSP detection\n\n"
   "-box : draws the CC search window on <prefix>_cc.ppm\n\n"
   "-W: outputs <prefix>_cc_witelson.ppm as well as 7 Witelson sub-areas\n\n"
   "-Hampel: outputs <prefix>_cc_hampel.ppm as well as 5 Hampel sub-areas\n\n"
   "-border : outputs <prefix>_cc_border.ppm\n\n"
   "-h : prints help message\n\n"
   "-n <# atlases> : number of atlases used (default: 49)\n\n"
   "-a <atlas> : atlas used (default: $ARTHOME/default_yuki_atlas.nii)\n\n"
   "-o <output prefix> : prefix for naming output files (default: same as input prefix)\n\n"
   "-csv <csvfile> : CC measurments (area, perimeter, etc.) will be appended to this file in\n"
   "comma-separated values (CSV) format (default: <output prefix>.csv)\n\n"
   "-i <subject volume> : the volume (short int NIFTI format) on which the corpus callosum\n"
   "is to be located\n\n");

   return;
}

//////////////////////////////////////////////////////////////////////////////////
int vertical_runs(short *cc, int nx, int ny, int i0)
{
   int sign_change=0;
   int runs=0;

   for(int j=0; j<ny-1; j++)
   {
      if( cc[j*nx + i0] != cc[(j+1)*nx + i0] ) sign_change++;
   }

   runs = sign_change/2;

   return( runs );
}

//////////////////////////////////////////////////////////////////////////////////
void cw_eight_neighbor(short *im, int nx, int ny, int del_i, int del_j, int &bi, int &bj, int &ci, int &cj)
{
   if( im[ (cj+del_j)*nx + (ci+del_i) ] == 1 )
   {
      ci += del_i;
      cj += del_j;
      //printf("ci=%d cj=%d\n", ci, cj);
      return;
   }
   else
   {
      bi = ci + del_i;
      bj = cj + del_j;
      // printf("bi=%d bj=%d\n", bi, bj);
   }


   if( del_i==1 && del_j==0)
   {
      cw_eight_neighbor(im, nx, ny, 1, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==1 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, 0, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==0 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, -1, 1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==1)
   {
      cw_eight_neighbor(im, nx, ny, -1, 0, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==0)
   {
      cw_eight_neighbor(im, nx, ny, -1, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==-1 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 0, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==0 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 1, -1, bi, bj, ci, cj);
      return;
   }
   else if( del_i==1 && del_j==-1)
   {
      cw_eight_neighbor(im, nx, ny, 1, 0, bi, bj, ci, cj);
      return;
   }

   return;
}

//////////////////////////////////////////////////////////////////////////////////
// NOTE: array t is expected to be a least size (n+1)
float compute_perimeter(float *x, float *y, float *t, int n)
{
   float delx, dely;

   t[0] = 0.0;

   for(int i=1; i<n; i++)
   {
      delx = x[i]-x[i-1];
      dely = y[i]-y[i-1];

      t[i] = t[i-1] + sqrtf( delx*delx + dely*dely );
   }

   // closing loop
   delx = x[n-1]-x[0];
   dely = y[n-1]-y[0];
   t[n] = t[n-1] + sqrtf( delx*delx + dely*dely );

   return(t[n]);
}
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// IMPORTANT: The posterior_point of the CC must be found correctly before running
// this function.
//////////////////////////////////////////////////////////////////////////////////
// nb: number of border pixels
float estimate_perimeter(short *msk, int nx, int ny, float dx, float dy, int *&cci, int*&ccj, int &nb)
{
   int ncc; // number of pixels comprising the CC
   int np;
   float *ccx, *ccy, *cct;
   int bi=0;
   int bj=0;
   int del_i, del_j;
   float perimeter;

   np = nx*ny;

   ncc=0;
   for(int i=0; i<np; i++)
   {
      if(msk[i]>0)
      {
         ncc++;
      }
   }

   // more memory than need because nb<ncc
   cci = (int *)calloc(ncc, sizeof(int));
   if(cci == NULL)
   {
      memory_allocation_error("cci");
   }

   ccj = (int *)calloc(ncc, sizeof(int));
   if(ccj == NULL)
   {
      memory_allocation_error("ccj");
   }

   ccx = (float *)calloc(ncc, sizeof(float));
   if(ccx == NULL)
   {
      memory_allocation_error("ccx");
   }

   ccy = (float *)calloc(ncc, sizeof(float));
   if(ccy == NULL)
   {
      memory_allocation_error("ccy");
   }

   cct = (float *)calloc(ncc, sizeof(float));
   if(cct == NULL)
   {
      memory_allocation_error("cct");
   }

   cci[0] = posterior_point[0];
   ccj[0] = posterior_point[1];

   // this is a critical choice.  It happens to work for the way the posterior point is found.  
   // It must be altered for a different starting point.
   del_i = 1; 
   del_j = 0;

   nb = 0;
   do
   {
      nb++;

      cw_eight_neighbor(msk, nx, ny, del_i, del_j, bi, bj, posterior_point[0], posterior_point[1]);

      del_i = bi - posterior_point[0];
      del_j = bj - posterior_point[1];

      cci[nb]=posterior_point[0];
      ccj[nb]=posterior_point[1];

   }
   while( posterior_point[0]!=cci[0] || posterior_point[1]!=ccj[0] );  // see if it has come full circle

   // convert from (i,j) to (x,y) coordinates
   for(int b=0; b<nb; b++)
   {
      ccx[b] = (cci[b]-(nx-1.0)/2.0)*dx;
      ccy[b] = (ccj[b]-(ny-1.0)/2.0)*dy;
   }

/*
This was not necessary and causing problems
   for(int b=0; b<nb; b++)
   {
      msk[ccj[b]*nx + cci[b]]=2;
   }
*/

   perimeter=compute_perimeter(ccx, ccy, cct, nb);

   free(ccx);
   free(ccy);
   free(cct);

   return( perimeter );
}

//////////////////////////////////////////////////////////////////////////////////

void estimate_witelson(short *msk, int nx, int ny, float dx, float dy, float *W)
{
   int np;
   int v;

   np = nx*ny;

   for(int i=0; i<=7; i++)
   {
      W[i]=0;
   }

   for(int i=0; i<nx; i++)
   for(int j=0; j<ny; j++)
   {
      v = j*nx + i;
      if(msk[v]>0 && i>=(posterior_point[0]-cc_length_fifth) ) // W7
      {
         W[7] += dx*dy;
      }
      else if(msk[v]>0 && i>=(posterior_point[0]-cc_length_third) ) // W6
      {
         W[6] += dx*dy;
      }
      else if(msk[v]>0 && i>=(posterior_point[0]-cc_length_half) ) // W5
      {
         W[5] += dx*dy;
      }
      else if(msk[v]>0 && i>=(anterior_point[0]+cc_length_third) ) // W4
      {
         W[4] += dx*dy;
      }
      else if(msk[v]>0 && i<=posterior_genu[0] ) // W2
      {
         W[2] += dx*dy;
      }
      else if(msk[v]>0 && j>=posterior_genu[1] ) // W1
      {
         W[1] += dx*dy;
      }
      else if(msk[v]>0) // W3
      {
         W[3] += dx*dy;
      }
   }
}

void estimate_hampel(short *msk, int nx, int ny, float dx, float dy, float *H)
{
   float d;
   float pi;
   float i0, j0;
   float theta;
   float costheta=0.0;

   pi = 4.0*atanf(1.0);

   i0 = hampel_origin[0];
   j0 = hampel_origin[1];

   int np;
   int v;

   np = nx*ny;

   for(int i=0; i<=5; i++)
   {
      H[i]=0;
   }

   for(int i=0; i<nx; i++)
   for(int j=0; j<ny; j++)
   {
      v = j*nx + i;

      d = sqrtf ( (i-i0)*(i-i0) + (j-j0)*(j-j0) );

      costheta = (i-i0)*hampel_axis[0]/d + (j-j0)*hampel_axis[1]/d ;
      if(costheta > 1.0 ) costheta=1.0;
      if(costheta < -1.0 ) costheta=-1.0;
      theta = acosf( costheta );

      if(msk[v]>0 && theta <= pi/5.0) // H5
      {
         H[5] += dx*dy;
      }
      else if(msk[v]>0 && theta <= 2.0*pi/5.0 ) // H4
      {
         H[4] += dx*dy;
      }
      else if(msk[v]>0 && theta <= 3.0*pi/5.0 ) // H3
      {
         H[3] += dx*dy;
      }
      else if(msk[v]>0 && theta <= 4.0*pi/5.0 ) // H2
      {
         H[2] += dx*dy;
      }
      else if(msk[v]>0 && theta <= pi ) // H1
      {
         H[1] += dx*dy;
      }
      else if(msk[v]>0) // H0
      {
         H[0] += dx*dy;
      }
   }
}
//////////////////////////////////////////////////////////////////////////////////

float estimate_area(short *msk, int nx, int ny, float dx, float dy)
{
   int ncc;
   int np;

   np = nx*ny;

   ncc=0;
   for(int i=0; i<np; i++)
   {
      if(msk[i]>0)
      {
         ncc++;
      }
   }

   return( ncc*dx*dy );
}
//////////////////////////////////////////////////////////////////////////////////

float compute_circularity(float area, float perimeter)
{
   float pi;
   float circularity=0.0;

   pi = 4.0*atanf(1.0);

   if(perimeter > 0.0)
   {
      circularity = 4.0*pi*area/(perimeter*perimeter);
   }

   return(circularity);
}

//////////////////////////////////////////////////////////////////////////////////

void find_posterior_genu(short *cc_est)
{
   int j0=0, j1=0;
   posterior_genu[0]=anterior_point[0];
   posterior_genu[1]=anterior_point[1];

   // starting a AC and going towards to posterior directions (i.e., decreasing i)
   for(int i=anterior_point[0]; i<=posterior_point[0]; i++)
   {
      // find the first i-position where there are two vertical runs
      if( vertical_runs(cc_est, NX, NY, i) == 2 )
      {
         posterior_genu[0] = i-1;

         for(int j=0; j<NY-1; j++)
         {
            if(cc_est[j*NX + i]>0 && cc_est[(j+1)*NX + i]==0) 
            {
               j0 = j;
               break;
            }
         }

         for(int j=NY-1; j>0; j--)
         {
            if(cc_est[j*NX + i]>0 && cc_est[(j-1)*NX + i]==0) 
            {
               j1 = j;
               break;
            }
         }

         break;
      }
   }

   posterior_genu[1] = (int)((j0+j1)/2.0);
}

//////////////////////////////////////////////////////////////////////////////////
void find_rostrum(short *cc_est)
{
   for(int i=NX/2; i>=0; i--)
   {

      // find the first i-position where there are two vertical runs
      if( vertical_runs(cc_est, NX, NY, i) == 2 )
      {
         int sign_changes=0;
         int j0=0, j1=0;
         for( int j=NY-1; j>0; j-- )
         {
            if( cc_est[j*NX + i] != cc_est[(j-1)*NX + i] )
            {
               sign_changes++;

               if(sign_changes == 1) 
               { 
                  j1 = j;
               }
               else if(sign_changes == 2) 
               { 
                  j0 = j-1;
                  break;
               }
            }
         }

         rostrum[0] = i;
         rostrum[1] = j1 - (j1-j0)/2;
         break;
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////
void find_inferior_genu(short *cc_est)
{
   int i0=0, j0=0, v=0;

   for(int i=0; i<NX/2; i++)
   for(int j=0; j<NY; j++)
   {
      v = j*NX + i;

      if( cc_est[v] > 0 )
      {
         if(j>j0) { j0=j; i0=i;}
      } 
   }

   inferior_genu[0] = i0;
   inferior_genu[1] = j0;
}

void find_inferior_splenium(short *cc_est)
{
   int i0=0, j0=0, v=0;

   for(int i=NX/2; i<NX; i++)
   for(int j=0; j<NY; j++)
   {
      v = j*NX + i;

      if( cc_est[v] > 0 )
      {
         if(j>j0) { j0=j; i0=i;}
      } 
   }

   /////////////////////////////////////////////////////////////////////
   // Ensure if there is a run of CC pixels at the bottom of the 
   // splenium the middle pixel is chosen.
   /////////////////////////////////////////////////////////////////////
   int i1, i2;

   i1=i2=i0;

   for(int i=NX/2; i<NX-1; i++)
   {
      v = j0*NX + i;

      if( cc_est[v] > 0 && cc_est[v-1] == 0) 
      {
         i1 = i;
         break;
      }
   }

   for(int i=NX-1; i>NX/2; i--)
   {
      v = j0*NX + i;

      if( cc_est[v] > 0 && cc_est[v+1] == 0) 
      {
         i2 = i;
         break;
      }
   }

   i0 = (i1 + i2)/2;
   /////////////////////////////////////////////////////////////////////
   
   inferior_splenium[0] = i0;
   inferior_splenium[1] = j0;
}

void find_hampel_coordinates()
{
   float ui=0.0, uj=0.0, d=0.0;
   float p0i=0.0, p0j=0.0; 
   float p1i=0.0, p1j=0.0; 
   float ri=0.0, rj=0.0;

   p0i = inferior_genu[0];
   p0j = inferior_genu[1];
   p1i = inferior_splenium[0];
   p1j = inferior_splenium[1];

   d = sqrtf( (p0i-p1i)*(p0i-p1i) + (p0j-p1j)*(p0j-p1j) );

   if(d!=0.0)
   {
      ui = (p1i-p0i)/d;
      uj = (p1j-p0j)/d;
   }

   ri = posterior_point[0]-p1i;
   rj = posterior_point[1]-p1j;
   p1i += (ri*ui + rj*uj)*ui;
   p1j += (ri*ui + rj*uj)*uj;

   ri = anterior_point[0]-p0i;
   rj = anterior_point[1]-p0j;
   p0i += (ri*ui + rj*uj)*ui;
   p0j += (ri*ui + rj*uj)*uj;

   hampel_origin[0] = p0i + (p1i-p0i)/2.0;
   hampel_origin[1] = p0j + (p1j-p0j)/2.0;

   hampel_axis[0] = ui;
   hampel_axis[1] = uj;
}

//////////////////////////////////////////////////////////////////////////////////
// This function splits the CC boundary (cci, ccj) into two parts (ubi, ubj) and
// (lbi, lbj) representing the upper and lower boundaries if the CC.  
// The first point of the low boundary coincides with the last point of the upper boundary.
// The first point of the upper boundary coindicides with the last point of the low boudary.
//////////////////////////////////////////////////////////////////////////////////
void separate_upper_and_lower_boundaries(int *&ubi, int *&ubj, int *&lbi, int *&lbj, int *cci, int *ccj, int nb)
{
   int b0=0, b1=0;

   // b0 is the index of the boundry point corresponding to the inferior splenium 
   // b1 is the index of the boundry point corresponding to the rostrum 
   for(int b=0; b<nb; b++)
   {
      if(ccj[b]==rostrum[1] && cci[b]==rostrum[0])
         b1=b;

      if(ccj[b]==inferior_splenium[1] && cci[b]==inferior_splenium[0])
         b0=b;
   }

   // nlb: number of points on the lower boundary of the CC defined as going
   // from b0 to b1 clockwise
   nlb = b1-b0+1;

   // nub: number of points on the upper boundary of the CC defined as going
   // from b1 to b0 clockwise
   nub = (b0 + 1) + (nb - b1);

   // (ubi, ubj) : coordinates of points on the upper boundary (0, 1, ..., nub-1)
   ubi = (int *)calloc(nub, sizeof(int));
   ubj = (int *)calloc(nub, sizeof(int));

   // (lbi, lbj) : coordinates of points on the lower boundary (0, 1, ..., nlb-1)
   lbi = (int *)calloc(nlb, sizeof(int));
   lbj = (int *)calloc(nlb, sizeof(int));

   // fill in (lbi, lbj)_b (b=0, ..., nlb-1) 
   // from (cci, ccj)_b0 to (cci, ccj)_b1
   for(int b=0; b<nlb; b++)
   {
      lbi[b] = cci[b+b0];
      lbj[b] = ccj[b+b0];
   }

   // fill in (lbi, lbj)_b (b=0, ..., nb-b1-1) 
   // from (cci, ccj)_b1 to (cci, ccj)_(nb-1)
   for(int b=0; b<(nb-b1); b++)
   {
      ubi[b] = cci[b+b1];
      ubj[b] = ccj[b+b1];
   }

   // fill in (lbi, lbj)_b (b=nb-b1, ..., nub-1) 
   // from (cci, ccj)_0 to (cci, ccj)_b0
   for(int b=(nb-b1); b<nub; b++)
   {
      ubi[b] = cci[b-nb+b1];
      ubj[b] = ccj[b-nb+b1];
   }

   return;
}
//////////////////////////////////////////////////////////////////////////////////
void find_medial_point(int i, int j, short *dist, short *cc)
{
   short mindist=32767;
   int ii=-1;
   int jj=-1;

   if(i<0 && j<0) return;

   medi[nm]=i;
   medj[nm]=j;
   nm++;

   cc[ j*NX + i ] = 2;

   if( i>0 && cc[ j*NX + i-1 ] == 1 )
   {
      cc[ j*NX + i-1 ] = 2;
      mindist = dist[ j*NX + i-1 ];
      ii=i-1; jj=j;
   }

   if( i>0 && j<NY-1 && cc[ (j+1)*NX + i-1 ] == 1 )
   {
      cc[ (j+1)*NX + i-1 ] = 2;
      if( dist[ (j+1)*NX + i-1 ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i-1 ];
         ii=i-1; jj=j+1;
      }
   }

   if( j<NY-1 && cc[ (j+1)*NX + i ] == 1 )
   {
      cc[ (j+1)*NX + i ] = 2;
      if( dist[ (j+1)*NX + i ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i ];
         ii=i; jj=j+1;
      }
   }

   if( i<NX-1 && j<NY-1 && cc[ (j+1)*NX + i+1 ] == 1 )
   {
      cc[ (j+1)*NX + i+1 ] = 2;
      if( dist[ (j+1)*NX + i+1 ] < mindist )
      {
         mindist = dist[ (j+1)*NX + i+1 ];
         ii=i+1; jj=j+1;
      }
   }

   if( i<NX-1 && cc[ j*NX + i+1 ] == 1 )
   {
      cc[ j*NX + i+1 ] = 2;
      if( dist[ j*NX + i+1 ] < mindist )
      {
         mindist = dist[ j*NX + i+1 ];
         ii=i+1; jj=j;
      }
   }

   if( i<NX-1 && j>0 && cc[ (j-1)*NX + i+1 ] == 1 )
   {
      cc[ (j-1)*NX + i+1 ] = 2;
      if( dist[ (j-1)*NX + i+1 ] < mindist )
      {
         mindist = dist[ (j-1)*NX + i+1 ];
         ii=i+1; jj=j-1;
      }
   }

   if( j>0 && cc[ (j-1)*NX + i ] == 1 )
   {
      cc[ (j-1)*NX + i ] = 2;
      if( dist[ (j-1)*NX + i ] < mindist )
      {
         mindist = dist[ (j-1)*NX + i ];
         ii=i; jj=j-1;
      }
   }

   if( i>0 && j>0 && cc[ (j-1)*NX + (i-1) ] == 1 )
   {
      cc[ (j-1)*NX + (i-1) ] = 2;
      if( dist[ (j-1)*NX + (i-1) ] < mindist )
      {
         mindist = dist[ (j-1)*NX + (i-1) ];
         ii=i-1; jj=j-1;
      }
   }

   find_medial_point(ii, jj, dist, cc);

   return;
}

void find_medial_axis(short *cc, short *dist)
{
   double dub, dlb;
   double dubmin, dlbmin;

   // (medi, medj) are coordinates of points on the medial axis
   medi= (int *)calloc(nub, sizeof(int));
   medj= (int *)calloc(nub, sizeof(int));

   // "distance" map pixel values equal the absolute value of the difference between
   // mininum distances of the pixel and the lower and upper CC boundaries
   dist = (short *)calloc(NX*NY, sizeof(short));

   // binarize cc (just in case)
   for(int v=0; v<NX*NY; v++) 
   if(cc[v]>0) { cc[v]=1; } else { cc[v]=0; }

   for(int i=0; i<NX; i++) 
   for(int j=0; j<NY; j++) 
   if( cc[j*NX + i] == 1)
   {
      dubmin = NX*dx+NY*dy; // can't be larger than this!
      for(int k=0; k<nub; k++)
      {
         dub = (ubi[k]-i)*dx*(ubi[k]-i)*dx + (ubj[k]-j)*dy*(ubj[k]-j)*dy;
         dub = sqrt(dub);
         if(dub<dubmin) dubmin=dub;
      }

      dlbmin = NX*dx+NY*dy; // can't be larger than this!
      for(int k=0; k<nlb; k++)
      {
         dlb = (lbi[k]-i)*dx*(lbi[k]-i)*dx + (lbj[k]-j)*dy*(lbj[k]-j)*dy;
         dlb = sqrt(dlb);
         if(dlb<dlbmin) dlbmin=dlb;
      }

      dist[j*NX + i]= (short)(100*fabs(dlbmin-dubmin));
   }
   else
   {
      dist[j*NX + i]= 0;
   }

   find_medial_point(inferior_splenium[0], inferior_splenium[1], dist, cc);
}

//////////////////////////////////////////////////////////////////////////////////
// nbp: number of boundary pixels. Equal nub, or nlb depending one whether we are
// interested in finding the upper or lower boundary crossing.
// (bi, bj) are boundary pixel locations
// (ui, uj) is unit normal poiting to the search direction.
// bp: index of the bp (in 1/2 increments)
// m: index of the medial axis pixel
void find_boundary_crossing(int i0, int j0, float &bp, short *cc, float ui, float uj, int *bi, int *bj, int nbp)
{
   float dpmax=-1.0;
   float dp;
   int indx;

   if(i0<1 || j0<1 || i0>=(NX-1) || j0>=(NY-1))
   {
      bp = -1.0;
      return;
   }

   if( cc[ j0*NX + i0 ] == 2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp; b++)
      {
         if(bi[b]==i0 && bj[b]==j0)
         {
            bp = b;
            break;
         }
      }
      return;
   }
   else if( cc[ (j0+1)*NX + (i0-1) ]==1 && cc[ j0*NX + (i0-1) ]==2 && cc[ (j0+1)*NX + i0 ]==2 )
   {
      bp = -1.0;

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0+1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0-1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0+1)*NX + (i0+1) ]==1 && cc[ j0*NX + (i0+1) ]==2 && cc[ (j0+1)*NX + i0 ]==2 )
   {
      bp = -1.0;

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0+1) && bj[b]==j0 )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0+1) )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0-1)*NX + (i0+1) ]==1 && cc[ (j0-1)*NX + i0 ]==2 && cc[ j0*NX + (i0+1) ]==2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0-1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0+1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }
   else if( cc[ (j0-1)*NX + (i0-1) ]==1 && cc[ (j0-1)*NX + i0 ]==2 && cc[ j0*NX + (i0-1) ]==2 )
   {
      bp = -1.0;
      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==i0 && bj[b]==(j0-1) )
         {
            bp = b;
            break;
         }
      }

      for(int b=0; b<nbp-1; b++)
      {
         if(bi[b]==(i0-1) && bj[b]==j0 )
         {
            bp += b;
            break;
         }
      }

      bp /= 2.0;
      return;
   }

   // NW pixel
   dp = ui*(-0.7071068) + uj*(-0.7071068);
   dpmax = dp; indx=0;

   // N pixel
   dp = -uj;
   if( dp > dpmax )
   {
      dpmax = dp; indx=1;
   }

   // NE pixel
   dp = ui*(0.7071068) + uj*(-0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=2;
   }

   // E pixel
   dp = ui;
   if( dp > dpmax )
   {
      dpmax = dp; indx=3;
   }

   // SE pixel
   dp = ui*(0.7071068) + uj*(0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=4;
   }

   // S pixel
   dp = uj;
   if( dp > dpmax )
   {
      dpmax = dp; indx=5;
   }

   // SW pixel
   dp = ui*(-0.7071068) + uj*(0.7071068);
   if( dp > dpmax )
   {
      dpmax = dp; indx=6;
   }

   // W pixel
   dp = -ui;
   if( dp > dpmax )
   {
      dpmax = dp; indx=7;
   }

   switch (indx) {
      case 0:
         find_boundary_crossing(i0-1, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 1:
         find_boundary_crossing(i0, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 2:
         find_boundary_crossing(i0+1, j0-1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 3:
         find_boundary_crossing(i0+1, j0, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 4:
         find_boundary_crossing(i0+1, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 5:
         find_boundary_crossing(i0, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 6:
         find_boundary_crossing(i0-1, j0+1, bp, cc, ui, uj, bi, bj, nbp);
         break;
      case 7:
         find_boundary_crossing(i0-1, j0, bp, cc, ui, uj, bi, bj, nbp);
         break;
   }

   return;
}

//////////////////////////////////////////////////////////////////////////////////

void find_thickness_profile(short *cc, const char *prefix)
{
   FILE *fp;
   char outputfile[1024]="";
   float dum=0.0;
   float *thickness_profile;
   short *arclength;
   float *tmp;

   ////////////////////////////////////////////////////////////////////
   // Find the normal vectors to the medial axis.
   // Program is designed such that the normal vectors always point 
   // towards the upper boundary of the CC.
   ////////////////////////////////////////////////////////////////////
   normi = (float *)calloc(nm, sizeof(float));
   normj = (float *)calloc(nm, sizeof(float));

   int m0, m1, dm=3;

   // find a normal vector for each pixel on the medial axis
   for(int m=0; m<nm; m++)
   {
      m0 = m - dm;
      m1 = m + dm;

      // take care of boundary conditions by reflection
      if( m0 < 0 ) m0 *= -1; 
      if( m1 >= nm ) m1 = 2*nm-2-m1;

      // at this point this is only the tangent line
      normi[m] = medi[m1] - medi[m0];
      normj[m] = medj[m1] - medj[m0];

      dum = sqrtf( normi[m]*normi[m] + normj[m]*normj[m] );

      if(dum>0.0)
      {
         normi[m] /= dum;
         normj[m] /= dum;
      }
      else
      {
         normi[m] = 0.0;
         normj[m] = 0.0;
      }

      // up to this point, normi and normj has been the tanget line
      // here, it converted to the normal line by a 90-deg rotation (i,j) --> (-j,i)
      dum = normi[m];
      normi[m] = -normj[m];
      normj[m] = dum;
   }
   ////////////////////////////////////////////////////////////////////

   thickness_profile = (float *)calloc(1001, sizeof(float));
   arclength = (short *)calloc(nm, sizeof(short));
   tmp = (float *)calloc(nm, sizeof(float));
   ////////////////////////////////////////////////////////////////////

   // indices of the nm boundary point where the nm normals to the
   // medial axis points cross the upper boundary
   // possible values are 0 to nub
   ubindx = (float *)calloc(nm, sizeof(float));

   // indices of the nm boundary point where the nm normals to the
   // medial axis points cross the lower boundary
   // possible values are 0 to nlb
   lbindx = (float *)calloc(nm, sizeof(float));

   // ensure that cc is 0 or 1
   for(int i=0; i<NX; i++)
   for(int j=0; j<NY; j++)
   if( cc[j*NX + i] > 0 )
   {
      cc[j*NX + i] = 1;
   }
   else
   {
      cc[j*NX + i] = 0;
   }

   // make the boundary points 2
   for(int b=0; b<nub; b++)
   {
      cc[ ubj[b]*NX + ubi[b] ] = 2;
   }

   for(int b=0; b<nlb; b++)
   {
      cc[ lbj[b]*NX + lbi[b] ] = 2;
   }

   for(int m=0; m<nm; m++)
   {
      find_boundary_crossing(medi[m], medj[m], ubindx[m], cc, normi[m], normj[m], ubi, ubj, nub);
      find_boundary_crossing(medi[m], medj[m], lbindx[m], cc, -normi[m], -normj[m], lbi, lbj, nlb);

      if( ubindx[m] < 0.0) ubindx[m]=0.0;
      if( lbindx[m] < 0.0) lbindx[m]=0.0;
   }

   float maxindx;
   float minindx;

   // ensure that ubindx values are in decreasing order
   // those points that violate this are set to zero
   maxindx=ubindx[0];
   for(int m=1; m<nm; m++)
   if( ubindx[m]>0.0 ) 
   {
      if( ubindx[m] <= maxindx )
      {
         maxindx = ubindx[m];
      }
      else
      {
         ubindx[m] = 0.0;
      }
   }

   // ensure that lbindx values are in increasing order
   // those points that violate this are set to zero
   minindx=lbindx[0];
   for(int m=1; m<nm; m++)
   if( lbindx[m]>0.0 ) 
   {
      if( lbindx[m] >= minindx )
      {
         minindx = lbindx[m];
      }
      else
      {
         lbindx[m] = 0.0;
      }
   }

   for(int m=1; m<(nm-1); m++)
   {
      if ( ubindx[m]==0.0 || lbindx[m]==0.0 )
      {
          ubindx[m]=lbindx[m]=0.0;
      }
   }

   // those points in ibindx and ubindx that were set to zero are filled
   // by interpolating between neighboring non-zero points
   // zero filing of lbindx and ubindx
   int count=0;
   for(int m=1; m<(nm-1); m++)
   {
      if ( ubindx[m]==0.0 && lbindx[m]==0.0 )
      {
         count=1;
         m0 = m;
         while( ubindx[m0 + count]==0.0 && lbindx[m0 + count]==0.0 )
         {
            count++;
         }
         m1 = m0+count;

         // interpolation
         for( int c=0; c<count; c++)
         {
            lbindx[m0+c] = lbindx[m0-1] + (lbindx[m1]-lbindx[m0-1])*(c+1)/(m1-m0+1);
            ubindx[m0+c] = ubindx[m0-1] + (ubindx[m1]-ubindx[m0-1])*(c+1)/(m1-m0+1);
         }

         // jump over the 00000 segment that was just processed
         m += count;
      }
   }

   tmp[0]=0.0;
   for(int m=1; m<nm; m++)
   {
      tmp[m] = tmp[m-1] + 
      sqrtf( (medi[m]-medi[m-1])*(medi[m]-medi[m-1])*dx*dx + (medj[m]-medj[m-1])*(medj[m]-medj[m-1])*dy*dy );
   }

   if( tmp[nm-1] > 0.0 )
   for(int m=0; m<nm; m++)
   {
      tmp[m] /= tmp[nm-1];
      arclength[m] = (short)(tmp[m]*10000.0 + 0.5);
   }

   float res;
   float lx, ly, ux, uy;
   for(int m=0; m<nm; m++)
   {
      res = lbindx[m] - (int)lbindx[m];

      if( (int)lbindx[m]< (nlb-1) )
      {
         lx = dx* ( lbi[ (int)lbindx[m] ]*(1.0-res) + lbi[ (int)lbindx[m] + 1]*res );
         ly = dy* ( lbj[ (int)lbindx[m] ]*(1.0-res) + lbj[ (int)lbindx[m] + 1]*res );
      }
      else
      {
         lx = dx*lbi[ (int)lbindx[m] ];
         ly = dy*lbj[ (int)lbindx[m] ];
      }

      res = ubindx[m] - (int)ubindx[m];

      if( (int)ubindx[m]< (nub-1) )
      {
         ux = dx * ( ubi[ (int)ubindx[m] ]*(1.0-res) + ubi[ (int)ubindx[m] + 1]*res );
         uy = dy * ( ubj[ (int)ubindx[m] ]*(1.0-res) + ubj[ (int)ubindx[m] + 1]*res );
      }
      else
      {
         ux = dx*ubi[ (int)ubindx[m] ];
         uy = dy*ubj[ (int)ubindx[m] ];
      }

      tmp[m] = sqrtf( (ux-lx)*(ux-lx) + (uy-ly)*(uy-ly) );
   }

   thickness_profile[0]=thickness_profile[1000]=0.0;
   for(int i=1; i<1000; i++)
   {
      for(int m=1; m<nm; m++)
      if( arclength[m] >= (i*10)  )
      {
         if( (arclength[m]- arclength[m-1]) > 0 )
            res = ((i*10.0) - arclength[m-1])/(arclength[m]- arclength[m-1]);
         else
            res = 0.0;
         thickness_profile[i] = tmp[m-1]*(1.0-res) + tmp[m]*res;
         break;
      }
   }

   free(tmp);

   tmp = smoothY(thickness_profile, 1, 1001, 1, 10.0);
   tmp[0]=tmp[1000]=0.0;

   sprintf(outputfile,"%s_TP.txt",prefix);
   fp = fopen(outputfile,"w");

   if(fp != NULL)
   {
      for(int i=0; i<1001; i+=10)
      {
         fprintf(fp,"%f\n",tmp[i]);
      }

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////////////
   free(arclength);
   free(tmp);
   free(normj);
   free(normi);
   free(thickness_profile);
   free(ubindx);
   free(lbindx);
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   FILE *fp;

   double startTime=0.0;
   double endTime=0.0;

   float max_t=0.0;

   short *atlas_cc=NULL;
   short *atlas_msp=NULL;

   int vox_offset=0;

   float W[8]; // areas of Witelson's subdivisions are stored in W[1], ... W[7]; W[0] is unused
   float H[6]; // areas of Hampel's subdivisions are stored in H[1], ... H[5]; H[0] is unused
   float CCarea;
   float CCperimeter;
   float CCcircularity;

   int number_of_atlases_used=49;
   float *corr=NULL;
   int *atlas_indx=NULL;
   int bbnx, bbny, bbnp;

   float *avg_warped_cc=NULL;
   float *dumf=NULL;
   short *cc_est=NULL;

   int count;
   short *atlas_msp_ptr;
   short *atlas_cc_ptr;

   char preselected_atlases[1024]="";
   char selected_atlases_file[1024]="";
   char subjectfile[1024]="";
   char csvfile[1024]="";

   int number_of_atlases_available;

   nifti_1_header atlas_hdr; // root only

   float sd;

   int N=11; // N=2*L+1

   char prefix[1024]="";
   char outputfile[1024]="";
   char inputfile[1024]="";
   char ccfile[1024]="";
   char msp_transformation_file[1024]="";

   short *trg=NULL;
   short *bbtrg;

   short *atlas=NULL;
   char atlasfile[1024]="yuki2.aux";

   int window_width=5;

   /////////////////////////////////////////////
   // MPI related
   int abort_flag=NO;
   int nproc; // Number of processes
   int root=0;
   int procID; // Process ID
   char procname[MPI_MAX_PROCESSOR_NAME];  // Processor name
   int procnamelen; // Length of the processor name

   /////////////////////////////////////////////

   // Initialize MPI
   MPI_Init(&argc, &argv);
   
   // Get process rank
   MPI_Comm_rank(MPI_COMM_WORLD, &procID);

   if(procID==root)
   {
      startTime = MPI_Wtime();
   }
   // Get total number of processes specificed at start of run
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);

   MPI_Get_processor_name(procname, &procnamelen);

   // by setting this option, the program will not output the *ACPC_axial.ppm and *ACPC_sagittal.ppm files
   // for the AC/PC detection program.
   opt_ppm=0;
   opt_txt=0;

   ////////////////////////////////////////////////////////////////////////
   // Only root process write messages
   if( argc==1 ) 
   {
      if( procID==root ) { print_help(); }
      MPI_Finalize(); // MPI needs to be always finalized before exit
      exit(0);
   }
   ////////////////////////////////////////////////////////////////////////

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 's':
            if( procID==root ) { print_secret_help(); }
            MPI_Finalize();
            exit(0);
         case 'V':
            if( procID==root ) { printf("Version 2.1 for Open MPI 1.6\n"); }
            MPI_Finalize();
            exit(0);
         case 'h':
            if( procID==root ) { print_help(); }
            MPI_Finalize();
            exit(0);
         case 'n':
            number_of_atlases_used=atoi(optarg);
            if(number_of_atlases_used<=0) number_of_atlases_used=49;
            break;
         case 'T':
            sprintf(msp_transformation_file,"%s",optarg);
            break;
         case 'a':
            sprintf(atlasfile,"%s",optarg);
            opt_a=YES;
            break;
         case 'L':
            sprintf(preselected_atlases,"%s",optarg);
            break;
         case 'c':
            sprintf(csvfile,"%s",optarg);
            break;
         case 'i':
            sprintf(subjectfile,"%s",optarg);
            break;
         case 'o':
            sprintf(prefix,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'm':
            opt_mrx=YES;
            break;
         case 'p':
            opt_ppm=YES;
            opt_txt=YES;
            break;
         case 'b':
            opt_box=YES;
            break;
         case 'W':
            opt_W=YES;
            break;
         case 'H':
            opt_H=YES;
            break;
         case 'B':
            opt_border=YES;
            break;
         case 'C':
            opt_cc=YES;
            sprintf(ccfile,"%s",optarg);
            break;
         case 't':
            max_t  = atof(optarg);
            break;
         case 'S':
            VSPS[0] = atoi(argv[optind-1]);
            VSPS[1] = atoi(argv[optind+0]);
            VSPS[2] = atoi(argv[optind+1]);
            opt_RP = NO;
            break;
         case 'A':
            AC[0] = atoi(argv[optind-1]);
            AC[1] = atoi(argv[optind+0]);
            AC[2] = atoi(argv[optind+1]);
            opt_AC= NO;
            break;
         case 'P':
            PC[0] = atoi(argv[optind-1]);
            PC[1] = atoi(argv[optind+0]);
            PC[2] = atoi(argv[optind+1]);
            opt_PC= NO;
            break;
         case '?':
            if( procID==root) { print_help(); }
            MPI_Finalize();
            exit(0);
      }
   }

   /////////////////////////////////////////////////////////////////////////////////////////////
   
   if(preselected_atlases[0]!='\0') 
   {
      int dum;
      fp = fopen(preselected_atlases,"r");
      if(fp==NULL) file_open_error(preselected_atlases);
      fscanf(fp,"%d\n",&number_of_atlases_used);
      for(int i=0; i<number_of_atlases_used; i++) fscanf(fp,"%d\n",&dum);
      fscanf(fp,"%f\n",&max_t);
      fclose(fp);
   }

   /////////////////////////////////////////////////////////////////////////////////////////////

   if(max_t<0.0 || max_t>100.0) max_t=0.0;

   /////////////////////////////////////////////////////////////////////////////////////////////
   // get the value of the ARTHOME environment variable
   // The getenv() function searches the environment list for a string that matches "ARTHOME".
   // It returns a pointer to the value in the environment, or NULL if there is no match.

   char *ARTHOME;  // full path of the directory of the ART software

   ARTHOME=getenv("ARTHOME");

   if(ARTHOME == NULL)
   {
      if(procID==root) { printf("The ARTHOME environment variable is not defined.\n"); }
      MPI_Finalize();
      exit(0);
   }

   if(opt_v && procID==root)
   {
      printf("ARTHOME = %s\n",ARTHOME);
   }
   /////////////////////////////////////////////////////////////////////////////////////////////
   
   ///////////////////////////////////////////////////////////////////////////////
   // read subject volume

   if(!opt_cc)
   {

      if( subjectfile[0]=='\0')
      {
         if(procID==root) { printf("Please specify a subject volume using -i argument.\n"); }
         MPI_Finalize();
         exit(0);
      }
      ///////////////////////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////////////////////
      // check to see if subjectfile appears to be a NIFTI image
      if( procID==root && not_magical_nifti(subjectfile) )
      {
         abort_flag = YES;
      }

      MPI_Bcast(&abort_flag, 1, MPI_INT, root, MPI_COMM_WORLD);
      if(abort_flag)
      {
         MPI_Finalize();
         exit(0);
      }
      ///////////////////////////////////////////////////////////////////////////////

      if( opt_v && procID==root)
      {
         printf("Subject volume = %s\n",subjectfile);
      }

      ///////////////////////////////////////////////////////////////////////////////
      // only the root process reads the subjectfile
   
      if( procID == root)
      {
         subject_volume=(short *)read_nifti_image(subjectfile, &sub_hdr);

         if(subject_volume == NULL)
         {
            printf("Error reading %s, aborting ...\n", subjectfile);
            abort_flag = YES;
         }

         if( !abort_flag )
         {
            Snx=sub_hdr.dim[1]; Sny=sub_hdr.dim[2]; Snz=sub_hdr.dim[3];
            Sdx=sub_hdr.pixdim[1]; Sdy=sub_hdr.pixdim[2]; Sdz=sub_hdr.pixdim[3];

            if(sub_hdr.datatype != DT_SIGNED_SHORT && sub_hdr.datatype != 512)
            {
               printf("\nSorry, this program currently only handles volume of datatype DT_SIGNED_SHORT (4)\n"
               "or DT_UINT16 (512). %s has datatype %d.\n.", subjectfile, sub_hdr.datatype);
               abort_flag = YES;
            }
         }

         if( !abort_flag && opt_v)
         {
            printf("Subject volume matrix size = %d x %d x %d (voxels)\n", Snx, Sny, Snz);
            printf("Subject voxel size = %8.6f x %8.6f x %8.6f (mm3)\n", Sdx, Sdy, Sdz);
         }
      }

      MPI_Bcast(&abort_flag, 1, MPI_INT, root, MPI_COMM_WORLD);
      if(abort_flag)
      {
         MPI_Finalize();
         exit(0);
      }
   }
   /////////////////////////////////////////////////////////////////////////////////////////

   if( prefix[0]=='\0')
   {
      size_t L;

      if(!opt_cc)
      {
         L = strlen( (const char *)subjectfile);

         if(L-4>0)
         {
            strncpy(prefix, (const char *)subjectfile, L-4);
            prefix[L-4]='\0';
         }
         else
         {
            if(procID==root) {printf("Please specify an output prefix using -o option.\n");}
            MPI_Finalize();
            exit(0);
         }
      }
      else if(opt_cc)
      {
         L = strlen( (const char *)ccfile);

         if(L-4>0)
         {
            strncpy(prefix, (const char *)ccfile, L-4);
            prefix[L-4]='\0';
         }
         else
         {
            if(procID==root) {printf("Please specify an output prefix using -o option.\n");}
            MPI_Finalize();
            exit(0);
         }
      }
   }

   if(procID==root && opt_v)
   {
      printf("Output prefix = %s\n",prefix);
   }

   /////////////////////////////////////////////////////////////////////////////////////////
   
   if(!opt_cc)
   {
      if(!opt_a)
      {
         sprintf(inputfile,"%s/%s",ARTHOME, atlasfile);
      }
      else
      {
         sprintf(inputfile,"%s", atlasfile);
      }

      if(procID==root && opt_v)
      {
         printf("Atlas file = %s\n", inputfile);
      }

      /////////////////////////////////////////////////////////////////////////////////////////

      // check to see if atlasfile appears to be a NIFTI image
      if( procID==root && not_magical_nifti(inputfile) )
      {
         abort_flag = YES;
      }

      MPI_Bcast(&abort_flag, 1, MPI_INT, root, MPI_COMM_WORLD);
      if(abort_flag)
      {
         MPI_Finalize();
         exit(0);
      }

      /////////////////////////////////////////////////////////////////////////////////////////
      // root process reads the atlas

      if( procID==root )
      {
         atlas=(short *)read_nifti_image(inputfile, &atlas_hdr);

         if(atlas==NULL)
         {
            printf("Error reading %s, aborting ...\n", inputfile);
            abort_flag = YES;
         }
      }
      MPI_Bcast(&abort_flag, 1, MPI_INT, root, MPI_COMM_WORLD);
      if(abort_flag)
      {
         MPI_Finalize();
         exit(0);
      }

      /////////////////////////////////////////////////////////////////////////////////////////

      if( procID == root )
      {
         bbnx=atlas_hdr.dim[1]; 
         bbny=atlas_hdr.dim[2];
         number_of_atlases_available = atlas_hdr.dim[3]/2;
         dx=atlas_hdr.pixdim[1]; 
         dy=atlas_hdr.pixdim[2]; 
         dz=atlas_hdr.pixdim[3];
         vox_offset = (int)atlas_hdr.vox_offset;
      }
      MPI_Bcast(&bbnx, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&bbny, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&dx, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
      MPI_Bcast(&dy, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
      MPI_Bcast(&dz, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
      MPI_Bcast(&number_of_atlases_available, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&vox_offset, 1, MPI_INT, root, MPI_COMM_WORLD);
      bbnp= bbnx*bbny;

      /////////////////////////////////////////////////////////////////////////////////////////

      if(procID==root && atlas_hdr.datatype != DT_SIGNED_SHORT && atlas_hdr.datatype != 512)
      {
         printf("\nSorry, this program currently only handles atlases of datatype DT_SIGNED_SHORT (4)\n"
         "or DT_UINT16 (512). %s has datatype %d.\n.", atlasfile, atlas_hdr.datatype);
         abort_flag = YES;
      }
      MPI_Bcast(&abort_flag, 1, MPI_INT, root, MPI_COMM_WORLD);
      if(abort_flag)
      {
         MPI_Finalize();
         exit(0);
      }

      /////////////////////////////////////////////////////////////////////////////////////////

      if(procID==root)
      {
         output_hdr = atlas_hdr;
         output_hdr.dim[3] = 1;
         output_hdr.dim[1] = NX;
         output_hdr.dim[2] = NY;
         sprintf(output_hdr.descrip,"Created by ART ccdetector");
      }

      /////////////////////////////////////////////////////////////////////////////////////////

      if(N<=0) N=11;
      if( (N%2)==0 ) N += 1;			// make sure it's odd
	
      if(window_width<=0) window_width=5;
      if( (window_width%2)==0 ) window_width+=1;			// make sure it's odd

      if( niter<=0 ) niter=4;
	
      // correlation window is (2*Lx+1)x(2*Ly+1)
      Lx = (N-1)/2;
      Ly = (int)(Lx*dx/dy + 0.5);
      N2 = (2*Lx+1)*(2*Ly+1);

      //Search window size is (2*Wx+1)*(2*Wy+1)
      Wx = (window_width-1)/2;
      Wy = (int)(Wx*dx/dy + 0.5);

      sd = 2.0*Wx;

      /////////////////////////////////////////////////////////////////////////////////////////
      // allocate memory for various arrays

      if(procID==root)
      {
         cc_est = (short *)calloc(NP, sizeof(short));
         dumf = (float *)calloc(bbnp, sizeof(float));
      }

      if(preselected_atlases[0]!='\0') 
      {
         if(procID==root)
            printf("Preselected atlases = %s\n", preselected_atlases);

         fp = fopen(preselected_atlases,"r");
         if(fp==NULL) file_open_error(preselected_atlases);
         fscanf(fp,"%d\n",&number_of_atlases_used);
         fclose(fp);
      }
   
      /////////////////////////////////////////////////////////////////////////////////////////
      // ensures that the number of atlases used does no exceed the number of atlases available
      
      if(number_of_atlases_available < number_of_atlases_used)
      {
         number_of_atlases_used = number_of_atlases_available;
      }

      if(procID==root && opt_v)
      {
         //printf("process %d sees:\n",procID);
         //printf("Atlas matrix size = %d x %d (pixels)\n", bbnx, bbny);
         //printf("Atlas voxel size = %8.6f x %8.6f x %8.6f (mm3)\n", dx, dy, dz);
         printf("Number of atlases available = %d\n", number_of_atlases_available);
         printf("Number of atlases used = %d\n", number_of_atlases_used);
      }

      /////////////////////////////////////////////////////////////////////////////////////////
      // All processes need to allocate memory for these

      Xwarp=(float *)calloc(bbnp,sizeof(float));
      Ywarp=(float *)calloc(bbnp,sizeof(float));
      Zwarp=(float *)calloc(bbnp,sizeof(float));
      ARobj=(float *)calloc(N2,sizeof(float));
      ARtrg=(float *)calloc(N2,sizeof(float));
      corr=(float *)calloc(number_of_atlases_available,sizeof(float));
      atlas_indx = (int *)calloc(number_of_atlases_available, sizeof(int));
      bbtrg = (short *)calloc(bbnp, sizeof(short));
      avg_warped_cc = (float *)calloc(bbnp, sizeof(float));
      atlas_cc = (short *)calloc(bbnp*number_of_atlases_used, sizeof(short));
      atlas_msp = (short *)calloc(bbnp*number_of_atlases_used, sizeof(short));

      /////////////////////////////////////////////////////////////////////////////////////////
      // trg image will be NX*NY (i.e., 512*512)
   
      if(procID==root)
      {
         if( msp_transformation_file[0]=='\0')
         {
            trg = find_subject_msp(subjectfile,prefix);
         }
         else
         {
            trg = find_subject_msp(subjectfile, prefix, msp_transformation_file);
         }
      }

      /////////////////////////////////////////////////////////////////////////////////////////
   
      if(procID==root)
      {
         count=0;
         for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
         for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
         {
            bbtrg[count] = trg[j*NX + i];
            count++;
         }

         for(int a=0; a<number_of_atlases_available; a++)
         {
            atlas_msp_ptr = atlas+(2*a)*bbnp;
            atlas_cc_ptr  = atlas+(2*a+1)*bbnp;

            atlas_indx[a] = a;
            corr[a] = pearsonCorrelation(bbtrg, atlas_msp_ptr, bbnp );

         }
         hpsort(number_of_atlases_available, corr, atlas_indx);

         // if atlases are preselected, read them into last elements of atlas_indx
         if(preselected_atlases[0]!='\0') 
         {
            fp = fopen(preselected_atlases,"r");
            if(fp==NULL) file_open_error(preselected_atlases);
            fscanf(fp,"%d\n",&number_of_atlases_used);
            for(int i=0; i<number_of_atlases_used; i++)
               fscanf(fp,"%d\n",&atlas_indx[number_of_atlases_available-1-i]);
            fclose(fp);
         }

         // saves the selected atlases
         sprintf(selected_atlases_file,"%s_A.txt",prefix);
         fp = fopen(selected_atlases_file, "w");
         if(fp==NULL) file_open_error(selected_atlases_file);
         fprintf(fp, "%d\n", number_of_atlases_used);
         for(int i=0; i<number_of_atlases_used; i++)
         {
            fprintf(fp, "%d\n", atlas_indx[number_of_atlases_available-1-i]);
         }
         fclose(fp);
      }

      if(procID==root)
      {
         int a;
         for(int i=0; i<number_of_atlases_used; i++)
         {
            a = atlas_indx[number_of_atlases_available-1-i];
            atlas_msp_ptr = atlas+(2*a)*bbnp;
            atlas_cc_ptr  = atlas+(2*a+1)*bbnp;

            for(int v=0; v<bbnp; v++)
            {
               atlas_msp[i*bbnp + v] = atlas_msp_ptr[v];
               atlas_cc[i*bbnp + v] = atlas_cc_ptr[v];
            }
         }
         free(atlas);
      }

      MPI_Bcast(bbtrg, bbnp, MPI_SHORT, root, MPI_COMM_WORLD);
      MPI_Bcast(atlas_indx, number_of_atlases_available, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(corr, number_of_atlases_available, MPI_FLOAT, root, MPI_COMM_WORLD);
      MPI_Bcast(atlas_msp, number_of_atlases_used*bbnp, MPI_SHORT, root, MPI_COMM_WORLD);
      MPI_Bcast(atlas_cc, number_of_atlases_used*bbnp, MPI_SHORT, root, MPI_COMM_WORLD);

      /////////////////////////////////////////////////////////////////////////////////////////
      //
//atlas_hdr.dim[3]=1;
//sprintf(outputfile,"ccw0.nii");
//save_nifti_image(outputfile, bbtrg, &atlas_hdr);
   
      //for(int i=1; i<number_of_atlases_used; i++)
      for(int i=0; i<number_of_atlases_used; i++)
      {
         int a;
         short *warped_cc;

         if(procID == i%nproc)
         {
            a = atlas_indx[number_of_atlases_available-1-i];
            atlas_msp_ptr = atlas_msp+i*bbnp;
            atlas_cc_ptr  = atlas_cc+i*bbnp;

//if(i>0 && i<6)
//{
//sprintf(outputfile,"ccw%d.nii",i);
//save_nifti_image(outputfile, atlas_msp_ptr, &atlas_hdr);
//sprintf(outputfile,"CCL%d.nii",i);
//save_nifti_image(outputfile, atlas_cc_ptr, &atlas_hdr);
//}

            //computeWarpField(atlas_msp, bbtrg, sd, bbnp, bbnx, bbny);
            computeWarpField(atlas_msp_ptr, bbtrg, sd, bbnp, bbnx, bbny);

            // ensure that cc is represented by a 0 or 100 binary image
            for(int v=0; v<bbnp; v++) 
            {
               if(atlas_cc_ptr[v] > 0)
               {
                  atlas_cc_ptr[v] = 100;
               }
               else
               {
                  atlas_cc_ptr[v] = 0;
               }
            }

//if(i>0 && i<6)
//{
//warped_cc=computeReslicedImage(atlas_msp_ptr, bbnx,bbny,1, dx,dy,dz, bbnx,bbny,1, dx,dy,dz, Xwarp, Ywarp, Zwarp);
//atlas_hdr.dim[3]=1;
//sprintf(outputfile,"warp%d.nii",i);
//save_nifti_image(outputfile, warped_cc, &atlas_hdr);
//free(warped_cc);
//}

            warped_cc=computeReslicedImage(atlas_cc_ptr, bbnx,bbny,1, dx,dy,dz, bbnx,bbny,1, dx,dy,dz, Xwarp, Ywarp, Zwarp);

            // copy warped cc back to its original memory
            for(int v=0; v<bbnp; v++) 
            {
               avg_warped_cc[v] += warped_cc[v];
            }

            free(warped_cc);

            if(opt_v)
            {
               printf("%d/%d: atlas #%03d corr=%6.4f processor=%d on %s\n", 
               i+1, number_of_atlases_used, a, corr[number_of_atlases_available-1-i], procID, procname);
            }
         }
      }

      ////////////////////////////////////////////////////////////

      if(procID==root)
      {
         MPI_Reduce(MPI_IN_PLACE, avg_warped_cc, bbnp, MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
      }
      else
      {
         MPI_Reduce(avg_warped_cc, avg_warped_cc, bbnp, MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
      }

      if(procID==root)
      {
         for(int v=0; v<bbnp; v++) 
         {
            avg_warped_cc[v] /= number_of_atlases_used;
         }
      }

//atlas_hdr.dim[3]=1;
//atlas_hdr.datatype=16; 
//sprintf(outputfile,"p.nii");
//save_nifti_image(outputfile, avg_warped_cc, &atlas_hdr);
//atlas_hdr.datatype=4; 
      ////////////////////////////////////////////////////////////

      if( procID==root )
      {
         float mean1, mean2;
         float ssd;
         int n1, n2;
         float *fdr; // Fisher's discriminant ratio

         copyarray(avg_warped_cc, dumf, bbnp);

         fdr = (float *)calloc(101, sizeof(float));

         for(int k=0; k<=100; k++)
         {
            float t;
            t = k*1.0;

            n1=n2=0;
            mean1=mean2=0.0;
            ssd=0.0;
            for(int i=0; i<bbnx; i++)
            {
               for(int j=0; j<bbny; j++)
               {
                  int voxel = j*bbnx + i;

                  if( avg_warped_cc[voxel] > t && dumf[voxel] > 0.0) 
                  { 
                     n1++;
                     mean1 += bbtrg[voxel];
                  }
                  else if ( dumf[voxel] > 0.0)
                  {
                     n2++;
                     mean2 += bbtrg[voxel];
                  }
               } // j
            } // i

            if( n1>0 ) mean1 /= n1; 
            if( n2>0 ) mean2 /= n2; 
            if( n1==0 || n2==0 ) mean1=mean2=0.0;

            for(int i=0; i<bbnx; i++)
            {
               for(int j=0; j<bbny; j++)
               {
                  int voxel = j*bbnx + i;

                  if( avg_warped_cc[voxel] > t && dumf[voxel] > 0.0) 
                  { 
                     ssd += (mean1-bbtrg[voxel])*(mean1-bbtrg[voxel]);
                  }
                  else if ( dumf[voxel] > 0.0)
                  {
                     ssd += (mean2-bbtrg[voxel])*(mean2-bbtrg[voxel]);
                  }
               } // j
            } // i

            if( n1==0 || n2==0 ) fdr[k]=0.0;

            if( (n1+n2-2.0) > 0 ) ssd /= (n1+n2-2.0);

            fdr[k] = sqrtf( (mean1 - mean2)*(mean1 - mean2) / ssd);

            // printf("%f %d %d %d mean1 = %f mean2 = %f fdr=%f\n", t, n1, n2, n1+n2, mean1, mean2, fdr[k]);
         }

         int maxidx=0;
         float maxfdr=0.0;
         for(int k=0; k<=100; k++)
         {
            if(fdr[k] > maxfdr) 
            {
               maxidx=k;
               maxfdr = fdr[k];
            }
         }
         free(fdr);

         //if(opt_v)
         //{
         //   printf("maxfdr=%f maxidx=%d\n",maxfdr, maxidx);
         //}

         if(max_t != 0.0)
         {
            if( opt_v )
            {
               printf("Manually selected threshold = %3.1f\n",max_t);
            }
         }
         else
         {
            max_t = maxidx*1.0;
            if( opt_v )
            {
               printf("Automatically selected threshold = %3.1f\n",max_t);
            }
         }

         // save the threshold used for label fusion
         fp = fopen(selected_atlases_file, "a");
         if(fp==NULL) file_open_error(selected_atlases_file);
         fprintf(fp, "%f\n", max_t);
         fclose(fp);

         for(int i=0; i<bbnx; i++)
         {
            for(int j=0; j<bbny; j++)
            {
               int voxel = j*bbnx + i;

               if( avg_warped_cc[voxel] > maxidx*1.0) 
               { 
                  cc_est[ (j+UPPER_LEFT_j)*NX + (i+UPPER_LEFT_i)]=1;
               }
               else 
               {
                  cc_est[ (j+UPPER_LEFT_j)*NX + (i+UPPER_LEFT_i)]=0;
               }
            } // j
         } // i
      }

   }

   if(opt_cc)
   {
      cc_est=(short *)read_nifti_image(ccfile, &output_hdr);
      dx=output_hdr.pixdim[1]; 
      dy=output_hdr.pixdim[1]; 
   }

   ////////////////////////////////////////////////////////////
   if(procID==root)
   {
      // locates the most anterior point of the CC
      for(int i=UPPER_LEFT_i; i<=LOWER_RIGHT_i; i++)
      {
         if( vertical_runs(cc_est, NX, NY, i) )
         {
            int j0, j1; 

            j1 = 0; j0 = NY;
            for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
            {
               if(cc_est[j*NX + i])
               {
                  if( j>j1 ) j1=j;
                  if( j<j0 ) j0=j;
               }
            }

            anterior_point[0]=i;
            anterior_point[1]=(j1+j0)/2;
         
            break;
         }
      }

      // locates the most posterior point of the CC
      for(int i=LOWER_RIGHT_i; i>=UPPER_LEFT_i; i--)
      {
         if( vertical_runs(cc_est, NX, NY, i) )
         {
            int j0, j1; 

            j1 = 0; j0 = NY;
            for(int j=UPPER_LEFT_j; j<=LOWER_RIGHT_j; j++)
            {
               if(cc_est[j*NX + i])
               {
                  if( j>j1 ) j1=j;
                  if( j<j0 ) j0=j;
               }
            }

            posterior_point[0]=i;
            posterior_point[1]=(j1+j0)/2;
         
            break;
         }
      }

      // locates the most posterior point of the genu (not rostrum)
      find_posterior_genu(cc_est); // required for Witelson

      find_inferior_genu(cc_est); // required for Hampel
      find_inferior_splenium(cc_est); // required for Hampel
      find_hampel_coordinates();
      find_rostrum(cc_est);

      // fix a rare degenerate case
      if (cc_est[ posterior_point[1]*NX + posterior_point[0] ]==0)
      {
         cc_est[ posterior_point[1]*NX + posterior_point[0] ]=1;
      }

      cc_length = posterior_point[0]-anterior_point[0];
      cc_length_fifth = (int)nearbyint( cc_length/5.0 );
      cc_length_third = (int)nearbyint( cc_length/3.0 );
      cc_length_half = (int)nearbyint( cc_length/2.0 );

/*
      if(opt_v)
      {
         printf("Most anterior point: %d, %d\n", anterior_point[0], anterior_point[1]);
         printf("Most posterior point: %d, %d\n", posterior_point[0], posterior_point[1]);
         printf("CC length = %f mm\n", cc_length*dx);
         printf("CC length / 5 = %d pixels \n", cc_length_fifth);
      }
*/
   }

   //////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   if(!opt_cc)
   {
      if(procID==root)
      {
         // write program output
         output_ppm(trg, cc_est, (const char *)prefix);

         output_hdr.pixdim[4]=ACi;
         output_hdr.pixdim[5]=PCi;
         output_hdr.pixdim[6]=ACx;
         output_hdr.pixdim[7]=PCx;

         sprintf(outputfile,"%s_cc.nii",prefix);
         save_nifti_image(outputfile, cc_est, &output_hdr);

         update_qsform( (const char *)outputfile, Tacpc );
         //////////////////////////////////////////////////////////////////////////
      }
   }
   
   if(procID==root)
   {
      CCarea = estimate_area(cc_est, NX, NY, dx, dy);
      CCperimeter = estimate_perimeter(cc_est, NX, NY, dx, dy, cci, ccj, nb);
 
      // separates the CC boundary (cci, ccj) into an upper part (ubi, ubj)
      // and a lower part (lbi, lbj)
      separate_upper_and_lower_boundaries(ubi, ubj, lbi, lbj, cci, ccj, nb);

      // finds (medi, medj) representing a CC medial axis
      find_medial_axis(cc_est, dist);
      //uncomment to save "distance map"
      //sprintf(outputfile,"tt.nii");
      //save_nifti_image(outputfile, dist, &output_hdr);
      
      CCcircularity = compute_circularity(CCarea, CCperimeter);

      find_thickness_profile(cc_est, prefix);

      if(!opt_cc)
         output_bounding_box_ppm(trg, cc_est, (const char *)prefix);

      if(opt_W)
      {
         estimate_witelson(cc_est, NX, NY, dx, dy, W);
      }

      if(opt_H)
      {
         estimate_hampel(cc_est, NX, NY, dx, dy, H);
      }

      if(csvfile[0]=='\0')
      {
         sprintf(csvfile,"%s.csv",prefix);
      }

      if(opt_v)
      {
         printf("CC area = %6.2f mm^2\n",CCarea);
         printf("CC perimeter = %6.2f mm\n",CCperimeter);
         printf("CC circularity = %8.6f\n",CCcircularity);
         printf("CC length = %5.1f mm\n",cc_length*dx);
         if(opt_W)
         {
            printf("Witelson's subdivisoins:\n");
            printf("\tW1 = %6.2f mm^2\n",W[1]);
            printf("\tW2 = %6.2f mm^2\n",W[2]);
            printf("\tW3 = %6.2f mm^2\n",W[3]);
            printf("\tW4 = %6.2f mm^2\n",W[4]);
            printf("\tW5 = %6.2f mm^2\n",W[5]);
            printf("\tW6 = %6.2f mm^2\n",W[6]);
            printf("\tW7 = %6.2f mm^2\n",W[7]);
         }
         if(opt_H)
         {
            printf("Hampels's subdivisoins:\n");
            printf("\tC1 = %6.2f mm^2\n",H[1]);
            printf("\tC2 = %6.2f mm^2\n",H[2]);
            printf("\tC3 = %6.2f mm^2\n",H[3]);
            printf("\tC4 = %6.2f mm^2\n",H[4]);
            printf("\tC5 = %6.2f mm^2\n",H[5]);
         }
      }

      if( csvfile[0]!='\0')
      {
         FILE *fp;

         if (checkFileExistence(csvfile)==0)
         {
            fp = fopen(csvfile,"a");
            if(fp==NULL) file_open_error(csvfile);
            fprintf(fp,"ID, CC_area, CC_perimeter, CC_circularity, CC_length");
            if(opt_W) 
            {
               fprintf(fp,", W1, W2, W3, W4, W5, W6, W7");
            }
            if(opt_H) 
            {
               fprintf(fp,", C1, C2, C3, C4, C5");
            }
            fprintf(fp,"\n");
            fclose(fp);
         }

         fp = fopen(csvfile,"a");
         if(fp==NULL) file_open_error(csvfile);
         fprintf(fp,"%s, ",prefix);
         fprintf(fp,"%6.2f, ",CCarea);
         fprintf(fp,"%6.2f, ",CCperimeter);
         fprintf(fp,"%8.6f, ",CCcircularity);
         fprintf(fp,"%5.1f",cc_length*dx);
         if(opt_W)
         {
            fprintf(fp,", %6.2f, ",W[1]);
            fprintf(fp,"%6.2f, ",W[2]);
            fprintf(fp,"%6.2f, ",W[3]);
            fprintf(fp,"%6.2f, ",W[4]);
            fprintf(fp,"%6.2f, ",W[5]);
            fprintf(fp,"%6.2f, ",W[6]);
            fprintf(fp,"%6.2f",W[7]);
         }
         if(opt_H)
         {
            fprintf(fp,", %6.2f, ",H[1]);
            fprintf(fp,"%6.2f, ",H[2]);
            fprintf(fp,"%6.2f, ",H[3]);
            fprintf(fp,"%6.2f, ",H[4]);
            fprintf(fp,"%6.2f",H[5]);
         }
         fprintf(fp,"\n");
         fclose(fp);
      }
   }

   ////////////////////////////////////////////////////////////

   endTime = MPI_Wtime();

   if(!opt_cc)
   {
      if(opt_v && procID==root)
      {
         printf("Processing time = %6.2lf\n",endTime-startTime);
      }
   }

   ////////////////////////////////////////////////////////////
   // free all allocated memory
   if(procID==root)
   {
      free(cc_est);
      free(cci);
      free(ccj);
   }

   if(!opt_cc)
   {
      if(procID==root)
      {
         free(dumf);
      }

      free(avg_warped_cc);
      free(Xwarp);
      free(Ywarp);
      free(Zwarp);
      free(ARobj);
      free(ARtrg);
      free(atlas_indx);
      free(corr);
      free(atlas_msp);
      free(atlas_cc);
   }

   ////////////////////////////////////////////////////////////
   // Terminate  MPI execution environment
   MPI_Finalize();
}
