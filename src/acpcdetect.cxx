// BUG FIX (2011-07-29): The default interpolation method was nearest neighbor. The -nn option would use
// trilinear interpolatin.  This was corrected.  Now the default interpolation method is trilinear and
// the -nn option forces the program to use nearest neighbor.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <volume.h>
#include <ctype.h>
#include <nifti1_io.h>
#include <niftiimage.h>
#include <babak_lib.h>

#define YES 1
#define NO 0

int opt;

static struct option options[] =
{
   {"-nn",0,'n'}, // does nearest neighbor interpolation

   {"-T", 1, 'r'},
   {"-M", 0, 'M'},  // make the midpoint between AC and PC the center of the FOV

   {"-sform",0,'s'},
   {"-qform",0,'q'},
   {"-noppm",0,'p'},
   {"-notxt",0,'t'},

   {"-AC",1,'A'},
   {"-PC",1,'P'},
   {"-VSPS",1,'S'},

   {"-V",0,'V'},
   {"-version",0,'V'},
   {"-Version",0,'V'},

   {"--T2",0,'T'},
   {"-T2",0,'T'},

   {"--rvsps",1,'0'},
   {"-rvsps",1,'0'},

   {"--rac",1,'1'},
   {"-rac",1,'1'},

   {"--rpc",1,'2'},
   {"-rpc",1,'2'},

   {"-D",0,'D'},

   {"-v",0,'v'},
   {"--verbose",0,'v'},
   {"-verbose",0,'v'},

   {"-m",1,'m'},
   {"--model",1,'m'},
   {"-model",1,'m'},

   {"-i",1,'i'},
   {"--image",1,'i'},
   {"-image",1,'i'},

   {"-orient",1,'O'},
   {"--orient",1,'O'},
   {"-O",1,'O'},

   {"-h",0,'h'},
   {"--help",0,'h'},
   {"-help",0,'h'},

   {"-o",1,'o'},
   {"-output",1,'o'},
   {"--output",1,'o'},

   {"-oo",1,'u'},

   {"-onx",1,'x'},
   {"-ony",1,'y'},
   {"-onz",1,'z'},

   {"-odx",1,'X'},
   {"-ody",1,'Y'},
   {"-odz",1,'Z'},

   {0,0,0}
};

int opt_M=NO;
int opt_nn=NO;

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
   printf("\nUsage: acpcdetect [-V/-version -h/-help -v/-verbose -m/-model <model> -rvsps <r> -rac <r> -rpc <r>]\n"
   "[-O/-orient <code> -D -M] [-o/-output <output volume> -oo <output orientation code>]\n"
   "[-onx <int> -ony <int> -onz <int> -odx <float> -ody <float> -odz <float> -sform -qform -noppm -notxt]\n"
   "[-AC <int> <int> <int>] [-T <filename>]\n"
   "-i/-image <input volume>\n\n"
   "Required arguments:\n"
   "-i or -image <input volume>: Input (test) image volume in NIFTI format of type `short'\n\n"
   "Optional arguments:\n"
   "-h or -help: Prints help information.\n"
   "-v or -verbose : Enables verbose mode\n"
   "-m or -model <model>: User-specified template model (default = $ARTHOME/T1acpc.mdl)\n"
   "-rvsps <r>: Search radius for VSPS (default = 50 mm)\n"
   "-rac <r>: Search radius for AC (default = 15 mm)\n"
   "-rpc <r>: Search radius for PC (default = 15 mm)\n"
   "-O or -orient <code>: Three-letter orientation code (See users' guide for more details). Examples:\n"
   "\t\tPIL for Posterior-Inferior-Left\n"
   "\t\tRAS for Right-Anterior-Superior\n"
   "-o or -output <filename>: If this option is present, the program outputs an AC/PC and MSP aligned image\n"
   "with the given filename.\n"
   "-oo <output orientation code>: Three-letter orientation code of the output image.  If this is not\n"
   "specified, then the output image will have the same orientation as the input image.\n"
   "-onx, -ony, -onz: x/y/z matrix dimensions of the output image (default=same as input image).\n"
   "-odx, -ody, -odz: x/y/z voxel dimensions of the output image (default=same as input image).\n"
   "-M: Make the midpoint between AC and PC the center of the output FOV.\n"
   "-T <filename>: Writes output rigid-body transformation matrix to specified <filename>.\n"
   "-D: Prints additional information\n\n"
   "Outputs:\n"
   "<input volume>_ACPC_sagittal.ppm: Sagittal view of the detected AC/PC locations\n"
   "<input volume>_ACPC_axial.ppm: Axial view of the detected AC/PC locations\n"
   "<input volume>_ACPC.txt: Stores the detected AC/PC coordinates and the estimated mid-sagittal plane\n\n"
   );

   exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void computeSiemensVisionOffsets(float *Tmsp, float *AC, float *PC)
{
	double pi;
   	double beta;
   	double alpha;
	float *invT;
	float ac[4], pc[4];
	float pcac[4];
	float dum;
	float n[4];
	float N[4];
	float TLAI[16];
	float Rx[16], Ry[16], R[16];

	pi=4.0*atan(1.0);

	invT=inv4(Tmsp);

	// assumes original is PIL, N[3]=0.0 is not a mistake, in fact very necessary
	N[0]=0.0; N[1]=0.0; N[2]=1.0; N[3]=0.0;
	multi(invT,4,4,N,4,1,N);
	normalizeVector(N,3);

	multi(invT, 4, 4,  AC, 4,  1, ac); // assumes original is PIL
	multi(invT, 4, 4,  PC, 4,  1, pc); // assumes original is PIL

	free(invT);

	pcac[0] = ac[0] - pc[0];
	pcac[1] = ac[1] - pc[1];
	pcac[2] = ac[2] - pc[2];

	normalizeVector(pcac,3);

	crossProduct(pcac,N,n);
	n[3]=1.0;

	inversePILtransform("LAI", TLAI);
	multi(TLAI, 4, 4, n, 4, 1, n);

   	beta = asin( -(float)n[0] );
   	alpha = asin( -n[1]/cos(beta) );

   	beta *= 180/pi;
   	alpha *= 180/pi;

	Rx[0]=1.0; Rx[1]=0.0; Rx[2]=0.0; Rx[3]=0.0;
	Rx[4]=0.0; Rx[5]=(float)cos(alpha*pi/180); Rx[6]=-(float)sin(alpha*pi/180); Rx[7]=0.0;
	Rx[8]=0.0; Rx[9]=(float)sin(alpha*pi/180); Rx[10]=(float)cos(alpha*pi/180); Rx[11]=0.0;
	Rx[12]=0.0; Rx[13]=0.0; Rx[14]=0.0; Rx[15]=1.0;

	Ry[0]=(float)cos(beta*pi/180); Ry[1]=0.0;	Ry[2]=-(float)sin(beta*pi/180); Ry[3]=0.0;
	Ry[4]=0.0; Ry[5]=1.0; Ry[6]=0.0; Ry[7]=0.0;
	Ry[8]=(float)sin(beta*pi/180); Ry[9]=0.0; Ry[10]=(float)cos(beta*pi/180); Ry[11]=0.0;
	Ry[12]=0.0; Ry[13]=0.0; Ry[14]=0.0; Ry[15]=1.0;

	multi(Rx,4,4,Ry,4,4,R);

	multi(TLAI, 4, 4,  ac, 4,  1, ac);
	multi(TLAI, 4, 4,  pc, 4,  1, pc);

	invT=inv4(R);
	multi(invT,4,4,ac,4,1,ac);
	multi(invT,4,4,pc,4,1,pc);
	free(invT);

	printf("\nSiemens Vision FOV offsets (assuming PIL orientation of structural scan):\n");
   	printf("Trans. to coronal = %lf deg.\nto sagittal = %lf deg.\n", alpha, beta);
	printf("Shift = %f mm\n",ac[2]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   int opt_D=NO;
   int opt_T2=NO;

   //////////////////////////////////////////////////
   float n[4];
   float d;
   //////////////////////////////////////////////////

   char transformation_filename[512];
   char outputfilename[512];
   char modelfile[1024];
   char imagefilename[512];
   char orientation[4];
   char output_orientation[4];

   float AC[4]={0.0, 0.0, 0.0, 1.0};
   float PC[4]={0.0, 0.0, 0.0, 1.0};
   float VSPS[4]={0.0, 0.0, 0.0, 1.0};
   float Tmsp[16]; // transforms volOrig to MSP aligned PIL orientation

   // searchradius[0] is for VSPS
   // searchradius[1] is for AC
   // searchradius[2] is for PC
   double searchradius[3]; // in units of mm

   int onx, ony, onz;
   float odx, ody, odz;

   nifti_1_header hdr; // 348 bytes
   ///////////////////////////////////////////////////////////////////////////////////////////////

   // It is very important to have these initializations
   onx=ony=onz=0;
   odx=ody=odz=0.0;
   outputfilename[0]='\0';
   transformation_filename[0]='\0';
   orientation[0]='\0';
   output_orientation[0]='\0';
   modelfile[0]='\0';
   imagefilename[0]='\0';
   searchradius[0] = 50.0;
   searchradius[1] = 15.0;
   searchradius[2] = 15.0;

   ///////////////////////////////////////////////////////////////////////////////////////////////

   while ((opt = getoption(argc, argv, options)) != -1 )
   {
      switch (opt)
      {
         case 'n':
            opt_nn=YES;
            break;
         case 'M':
            opt_M=YES;
            break;
         case 'A':
            AC[0] = atoi(argv[optind-1]);
            AC[1] = atoi(argv[optind+0]);
            AC[2] = atoi(argv[optind+1]);
            opt_AC = NO;
            break;
         case 'P':
            PC[0] = atoi(argv[optind-1]);
            PC[1] = atoi(argv[optind+0]);
            PC[2] = atoi(argv[optind+1]);
            opt_PC = NO;
            break;
         case 'S':
            VSPS[0] = atoi(argv[optind-1]);
            VSPS[1] = atoi(argv[optind+0]);
            VSPS[2] = atoi(argv[optind+1]);
            opt_RP = NO;
            break;
         case 't':
            opt_txt = NO;
            break;
         case 'p':
            opt_ppm = NO;
            break;
         case 's':
            opt_sform = YES;
            break;
         case 'q':
            opt_qform = YES;
            break;
         case 'V':
            printf("2011-04-05\n");
            exit(0);
         case 'x':
            onx = atoi(optarg);
            break;
         case 'y':
            ony = atoi(optarg);
            break;
         case 'z':
            onz = atoi(optarg);
            break;
         case 'X':
            odx = atof(optarg);
            break;
         case 'Y':
            ody = atof(optarg);
            break;
         case 'Z':
            odz = atof(optarg);
            break;
         case 'T':
			opt_T2=YES;
            break;
         case 'o':
            sprintf(outputfilename,"%s",optarg);
            break;
         case 'r':
            sprintf(transformation_filename,"%s",optarg);
            break;
         case 'i':
            sprintf(imagefilename,"%s",optarg);
            break;
         case '0':
            searchradius[0] = atof(optarg);
            break;
         case '1':
            searchradius[1] = atof(optarg);
            break;
         case '2':
            searchradius[2] = atof(optarg);
            break;
         case 'D':
            opt_D=YES;
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'O':
            sprintf(orientation,"%s",optarg);
            break;
         case 'u':
            sprintf(output_orientation,"%s",optarg);
            break;
         case 'm':
            sprintf(modelfile,"%s",optarg);
            break;
         case '?':
            print_help_and_exit();
      }
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////

   if( imagefilename[0]=='\0' )
   {
      printf("Please specify an input image using -i or -image argument.\n");
      exit(0);
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////

   if( not_magical_nifti(imagefilename) )
   {
      exit(0);
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("\n-------------------------------------------------------\n");
   }

   if(searchradius[0]<=0 || searchradius[0]>200.0) searchradius[0]=50.0;
   if(searchradius[1]<=0 || searchradius[1]>100.0) searchradius[1]=15.0;
   if(searchradius[2]<=0 || searchradius[2]>100.0) searchradius[2]=15.0;

   if(opt_D)
   {
      printf("\nVSPS search radius = %5.1lf mm\n",searchradius[0]);
      printf("AC search radius = %5.1lf mm\n",searchradius[1]);
      printf("PC search radius = %5.1lf mm\n",searchradius[2]);
   }

   detect_AC_PC_MSP(imagefilename, orientation, modelfile, searchradius, AC, PC, VSPS, Tmsp, opt_D, opt_v, opt_T2);

   if(opt_sform || opt_qform)
   {
      float T_ijk2xyz[16];
      float Tacpc[16];
      float PIL2RAS[16];
      float ac[4], pc[4];
      int N; // number of images at the end of the command line
      DIM input_dim;

      hdr = read_NIFTI_hdr(imagefilename);
      input_dim.nx = hdr.dim[1];
      input_dim.ny = hdr.dim[2];
      input_dim.nz = hdr.dim[3];
      input_dim.dx = hdr.pixdim[1];
      input_dim.dy = hdr.pixdim[2];
      input_dim.dz = hdr.pixdim[3];

      // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space and store
      // them in ac and pc variables.
      for(int i=0; i<4; i++) ac[i] = AC[i];
      for(int i=0; i<4; i++) pc[i] = PC[i];
      orig_ijk_to_pil_xyz(Tmsp, input_dim, ac, pc);

      // PIL2RAS takes points from PIL space to RAS space
      inversePILtransform("RAS", PIL2RAS);

      ACPCtransform(Tacpc, Tmsp, ac, pc, 1);

      // after this, Tacpc takes (x,y,z) points from original space to AC/PC-aligned RAS space
      multi(PIL2RAS, 4, 4,  Tacpc, 4,  4, Tacpc);

      ijk2xyz(T_ijk2xyz, input_dim.nx, input_dim.ny, input_dim.nz, input_dim.dx, input_dim.dy, input_dim.dz);

      // after this, Tacpc takes (i,j,k) points from original space to (x,y,z) points in AC/PC-aligned RAS space
      multi(Tacpc, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      // following lines were for testing only, but keep them for future
      // printMatrix(Tacpc, 4, 4, "sform", NULL);
      //multi(Tacpc, 4, 4,  AC, 4,  1, ac);
      //multi(Tacpc, 4, 4,  PC, 4,  1, pc);
      //printMatrix(ac, 4, 1, "ac", NULL);
      //printMatrix(pc, 4, 1, "pc", NULL);

      update_qsform( (const char *)imagefilename, Tacpc );

      N = argc-optind;

      for(int i=0; i<N; i++)
      {
         if(opt_v)
         {
            printf("\nCopying sform and/or qfrom information from %s to %s ...\n", imagefilename, (argv+optind)[i]);
         }

         if( !not_magical_nifti( (argv+optind)[i] ) )
         {
            update_qsform( (const char *)imagefilename, (const char*)(argv+optind)[i]);
         }
      }
   }

   if(transformation_filename[0]!='\0')
   {
      float Tacpc[16];
      float ac[4], pc[4];
      DIM input_dim;
      float PIL2OUT[16];
      float OUT2PIL[16];
      float PIL2RAS[16];
      float T_ijk2xyz[16];
      FILE *fp;

      hdr = read_NIFTI_hdr(imagefilename);
      input_dim.nx = hdr.dim[1];
      input_dim.ny = hdr.dim[2];
      input_dim.nz = hdr.dim[3];
      input_dim.dx = hdr.pixdim[1];
      input_dim.dy = hdr.pixdim[2];
      input_dim.dz = hdr.pixdim[3];

      // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
      for(int i=0; i<4; i++) ac[i] = AC[i];
      for(int i=0; i<4; i++) pc[i] = PC[i];
      orig_ijk_to_pil_xyz(Tmsp, input_dim, ac, pc);

      // PIL2OUT takes points from PIL space to output space
      if(output_orientation[0]!='\0' && isOrientationCodeValid(output_orientation) )
      {

         inversePILtransform(output_orientation, PIL2OUT);
      }
      else
      {

         inversePILtransform(orientation, PIL2OUT);
      }

      if(opt_M)
      {
         ACPCtransform(Tacpc, Tmsp, ac, pc, 0);
      }
      else
      {
         ACPCtransform(Tacpc, Tmsp, ac, pc, 1);
      }

      multi(PIL2OUT, 4, 4,  Tacpc, 4,  4, Tacpc);

      fp = fopen(transformation_filename,"w");
      if(fp==NULL) file_open_error(transformation_filename);
      printMatrix(Tacpc, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
      fclose(fp);
   }

   if(outputfilename[0]!='\0')
   {
      float Tacpc[16];
      DIM input_dim, output_dim;
      short *volOrig, *tmp;
      float PIL2OUT[16];
      float OUT2PIL[16];
      float PIL2RAS[16];
      float T_ijk2xyz[16];
      float *invT;

      volOrig = readNiftiImage(imagefilename, &input_dim, 0);

      if(onx<=0) onx=input_dim.nx;
      if(ony<=0) ony=input_dim.ny;
      if(onz<=0) onz=input_dim.nz;
      if(odx<=0.0) odx=input_dim.dx;
      if(ody<=0.0) ody=input_dim.dy;
      if(odz<=0.0) odz=input_dim.dz;

      output_dim.nx = onx;
      output_dim.ny = ony;
      output_dim.nz = onz;
      output_dim.dx = odx;
      output_dim.dy = ody;
      output_dim.dz = odz;

      // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
      orig_ijk_to_pil_xyz(Tmsp, input_dim, AC, PC);

      // PIL2OUT takes points from PIL space to output space
      if(output_orientation[0]!='\0' && isOrientationCodeValid(output_orientation) )
      {
         if(opt_v)
         {
            printf("\nOutput image orientation = %s\n", output_orientation);
         }
         inversePILtransform(output_orientation, PIL2OUT);
      }
      else
      {
         if(opt_v)
         {
            printf("\nOutput image orientation = %s\n", orientation);
         }
         inversePILtransform(orientation, PIL2OUT);
      }

      if(opt_M)
      {
         ACPCtransform(Tacpc, Tmsp, AC, PC, 0);
      }
      else
      {
         ACPCtransform(Tacpc, Tmsp, AC, PC, 1);
      }

      multi(PIL2OUT, 4, 4,  Tacpc, 4,  4, Tacpc);

      invT = inv4(Tacpc);
      if(!opt_nn)
      {
         tmp = resliceImage(volOrig,input_dim,output_dim,invT,LIN);
      }
      else
      {
         tmp = resliceImage(volOrig,input_dim,output_dim,invT,NEARN);
      }

      free(invT);

      hdr = read_NIFTI_hdr(imagefilename);
      hdr.pixdim[1]=odx; hdr.pixdim[2]=ody; hdr.pixdim[3]=odz;
      hdr.dim[1]=onx; hdr.dim[2]=ony; hdr.dim[3]=onz;
      hdr.magic[0]='n'; hdr.magic[1]='+'; hdr.magic[2]='1';
      sprintf(hdr.descrip,"Created by ART acpcdetect");
      save_nifti_image(outputfilename, tmp, &hdr);

      //////////////////////////////////////////////////////////////////////////////////
      // This part of the code adjusts the SFORM matrix of the output image

      if(output_orientation[0]!='\0' && isOrientationCodeValid(output_orientation) )
      {
         PILtransform(output_orientation, OUT2PIL);
      }
      else
      {
         PILtransform(orientation, OUT2PIL);
      }

      inversePILtransform("RAS", PIL2RAS);

      // NOTE: PIL orientation is a stepping stone here
      // Tacpc will transfrom for OUTPUT orientation to RAS orientation
      multi(PIL2RAS, 4, 4,  OUT2PIL, 4,  4, Tacpc);

      ijk2xyz(T_ijk2xyz, onx, ony, onz, odx, ody, odz);
      multi(Tacpc, 4, 4,  T_ijk2xyz, 4,  4, Tacpc);

      update_qsform( (const char *)outputfilename, Tacpc );
      //////////////////////////////////////////////////////////////////////////////////

      delete volOrig;
      delete tmp;
   }

   {
      compute_MSP_parameters_from_Tmsp(Tmsp, n, &d);
   }

   if(opt_v)
   {
      printf("\nEstimated mid-sagittal plane: (%8.7fx) + (%8.7fy) + (%8.7fz) = %8.5f (mm)\n", n[0],n[1],n[2],d);
   }

   if(opt_v)
   {
      printf("\nTest image: %s\n",imagefilename);
      printf("\nInput image orientation: %s\n",orientation);
      if( modelfile[0]!='\0' )
         printf("\nUsing model file: %s\n",modelfile);

      if(outputfilename[0]!='\0')
         printf("\nOutput image: %s\n",outputfilename);
   }
}
