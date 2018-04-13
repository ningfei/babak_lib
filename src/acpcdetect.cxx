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
   {"--standard", 0, 'S'},
   {"-standard", 0, 'S'},
   {"-s", 0, 'S'},
   {"-reorient", 1, 'R'},
   {"--reorient", 1, 'R'},
   {"-lm", 1, 'L'},
   {"-landmarks", 1, 'L'},
   {"--landmarks", 1, 'L'},
   {"-T", 1, 'r'},
   {"-centerAC", 0, 'M'},
   {"--centerAC", 0, 'M'},
//   {"-sform",0,'s'},
//   {"-qform",0,'q'},
   {"-noppm",0,'N'},
   {"-notxt",0,'t'},
   {"-V",0,'V'},
   {"-version",0,'V'},
   {"-Version",0,'V'},
//   {"--T2",0,'T'},
//   {"-T2",0,'T'},
   {"--rvsps",1,'0'},
   {"-rvsps",1,'0'},
   {"--rac",1,'1'},
   {"-rac",1,'1'},
   {"--rpc",1,'2'},
   {"-rpc",1,'2'},
   {"-v",0,'v'},
   {"--verbose",0,'v'},
   {"-verbose",0,'v'},
   {"-D",0,'v'},
//   {"-m",1,'m'},
//   {"--model",1,'m'},
//   {"-model",1,'m'},
   {"-i",1,'i'},
   {"--input",1,'i'},
   {"-input",1,'i'},
   {"-input_orient",1,'O'},
   {"--input_orient",1,'O'},
   {"-io",1,'O'},
   {"-h",0,'h'},
   {"--help",0,'h'},
   {"-help",0,'h'},
   {"-o",1,'o'},
   {"-output",1,'o'},
   {"--output",1,'o'},
   {"-oo",1,'u'},
   {"-output_orient",1,'u'},
   {"--output_orient",1,'u'},
   {"-onx",1,'x'},
   {"-ony",1,'y'},
   {"-onz",1,'z'},
   {"-nx",1,'x'},
   {"-ny",1,'y'},
   {"-nz",1,'z'},
   {"-odx",1,'X'},
   {"-ody",1,'Y'},
   {"-odz",1,'Z'},
   {"-dx",1,'X'},
   {"-dy",1,'Y'},
   {"-dz",1,'Z'},
   {0,0,0}
};

int opt_nn=NO;

//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit()
{
   printf("\nUsage: acpcdetect [-V/-version -h/-help -v/-verbose -rvsps <r> -rac <r> -rpc <r>]\n"
   "[--input_orient <code>] [-o/-output <output volume> -oo <output orientation code>]\n"
   "[-nx <int> -ny <int> -nz <int> -dx <float> -dy <float> -dz <float> -noppm -notxt]\n"
   "[-T <filename>] [-reorient Y/N (default=Y)]\n"
   "-i/-image <input volume>\n\n"
   "Required arguments:\n"
   "-i or -image <input volume>: Input (test) image volume in NIFTI format of type `short'\n\n"
   "Optional arguments:\n"
   "-h or -help: Prints help information.\n"
   "-v or -verbose : Enables verbose mode\n"
   "-noppm : Prevents outputting *.ppm images\n"
   "-notxt: Prevents outputting *_ACPC.txt\n"
//   "-m or -model <model>: User-specified template model (default = $ARTHOME/T1acpc.mdl)\n"
   "-rvsps <r>: Search radius for VSPS (default = 50 mm)\n"
   "-rac <r>: Search radius for AC (default = 15 mm)\n"
   "-rpc <r>: Search radius for PC (default = 15 mm)\n"
   "--input_orient <code>: Three-letter orientation code (See users' guide for more details). Examples:\n"
   "\t\tPIL for Posterior-Inferior-Left\n"
   "\t\tRAS for Right-Anterior-Superior\n"
   "-o or -output <filename>: If this option is present, the program outputs an AC/PC and MSP aligned image\n"
   "with the given filename.\n"
   "-oo <output orientation code>: Three-letter orientation code of the output image.  If this is not\n"
   "specified, then the output image will have the same orientation as the input image.\n"
   "-nx, -ny, -nz: x/y/z matrix dimensions of the output image (default=same as input image).\n"
   "-dx, -dy, -dz: x/y/z voxel dimensions of the output image (default=same as input image).\n"
   "-centerAC: Make AC the center of the output FOV.\n"
   "-T <filename>: Writes output rigid-body transformation matrix to specified <filename>.\n\n"
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
  char opt_standard=NO;
  char reorient='Y';
  float *invT;
  float Ttmp[16];
  FILE *fp;
  short *input_image;
  short *output_image;
  DIM input_dim, output_dim;
  nifti_1_header input_hdr, output_hdr; // 348 bytes
  float Tijk2xyz[16];
  float PIL2RAS[16];
  float PIL2OUT[16];
  float OUT2PIL[16];
//  int opt_T2=NO;
  float n[4]; 
  float d;
  char outputtransformationpath[1024]="";
  char landmarksfilepath[512]="";
//  char modelfile[1024]="";
  char inputimagepath[1024]="";
  char inputimagename[512]="";
  char inputimagedir[512]="";  // important to initialize to ""
  char outputimagepath[1024]="";
  char input_orient[4]="";
  char output_orient[4]="";
  float Tout[16]; // transforms input_image to the specified output orientation
  float TPIL[16]; // transforms input_image to PIL orientation

  // It is very important to have these initializations.
  int onx=0, ony=0, onz=0;
  float odx=0.0, ody=0.0, odz=0.0;

  // opt_CENTER_AC=NO means that by default the mid-point between AC and PC is set to the FOV
  // center.  If --centerAC is selected, then the AC is made the FOV center.
  // The opt_CENTER_AC variable is defined in PILtransform.cpp and made global in babak_lib.h
  opt_CENTER_AC=NO; 

  while ((opt = getoption(argc, argv, options)) != -1 )
  {
      switch (opt) 
      {
         case 'n':
            opt_nn=YES;
            break;
         case 'M':
            opt_CENTER_AC=YES;
            break;
         case 'S':
            opt_standard=YES;
            break;
         case 't':
            opt_txt = NO;
            break;
         case 'N':
            opt_ppm = NO;
            break;
//         case 's':
//            opt_sform = YES;
//            break;
//         case 'q':
//            opt_qform = YES;
//            break;
         case 'V':
            printf("April 10, 2018\n");
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
//         case 'T':
//            opt_T2=YES;
//            break;
         case 'R':
            reorient=optarg[0];
            break;
         case 'L':
            sprintf(landmarksfilepath,"%s",optarg);
            break;
         case 'o':
            sprintf(outputimagepath,"%s",optarg);
            break;
         case 'r':
            sprintf(outputtransformationpath,"%s",optarg);
            break;
         case 'i':
            sprintf(inputimagepath,"%s",optarg);
            break;
         case '0':
            searchradius[0] = atof(optarg); // searchradius[0] is for VSPS
            break;
         case '1':
            searchradius[1] = atof(optarg); // searchradius[1] is for AC
            break;
         case '2':
            searchradius[2] = atof(optarg); // searchradius[2] is for PC
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'h':
            print_help_and_exit();
            break;
         case 'O':
            sprintf(input_orient,"%s",optarg);
            break;
         case 'u':
            sprintf(output_orient,"%s",optarg);
            break;
//         case 'm':
//            sprintf(modelfile,"%s",optarg);
//            break;
         case '?':
            print_help_and_exit();
      }
  }

  if(searchradius[0]<=0 || searchradius[0]>200.0) searchradius[0]=50.0;
  if(searchradius[1]<=0 || searchradius[1]>100.0) searchradius[1]=15.0;
  if(searchradius[2]<=0 || searchradius[2]>100.0) searchradius[2]=15.0;

  if( inputimagepath[0]=='\0' )
  {
    printf("Please specify an input image using: -i <inputimage.nii>\n");
    exit(0);
  }
  if(opt_v) printf("Input image: %s\n",inputimagepath);

  if( outputimagepath[0]=='\0' )
  {
    printf("Please specify an output image using: -o <outputimage.nii>\n");
    exit(0);
  }
  if(opt_v) printf("Output image: %s\n",outputimagepath);

  // determine input image filename without the .nii suffix
  if( niftiFilename(inputimagename, inputimagepath)==0 ) { exit(1); }

  // determine input image directory
  getDirectoryName(inputimagepath, inputimagedir);

  if(outputtransformationpath[0]=='\0')
  {
    sprintf(outputtransformationpath,"%s/%s.mrx",inputimagedir,inputimagename);
  }
  if(opt_v) printf("Output transformation matrix: %s\n",outputtransformationpath);

  if(opt_v && landmarksfilepath[0]!='\0') 
  {
    printf("Manually specified landmarks: %s\n",landmarksfilepath);
  }

  input_image = (short  *)read_nifti_image(inputimagepath, &input_hdr);
  if(input_image==NULL)
  {
    printf("Error reading %s, aborting ...\n", inputimagepath);
    exit(1);
  }
  set_dim(input_dim, input_hdr);

  // if input orientation is specified using --input_orient option, make sure it's valid
  if(input_orient[0]!='\0' && isOrientationCodeValid(input_orient)==0)
  {
    printf("Error: %s is not a valid orientation code, aborting ...\n",input_orient);
    exit(0);
  }

  // if input orientation is not specified using --input_orient option, read it from input image
  if(input_orient[0]=='\0')
  {
    getNiftiImageOrientation(inputimagepath, input_orient);
    if(opt_v) printf("Input image orientation: %s\n",input_orient);
  } else {
    if(opt_v) printf("Input image orientation manually specified as: %s\n",input_orient);
  }

  // if ouput orientation is specified using -oo option, make sure it's valid
  if(output_orient[0]!='\0' && isOrientationCodeValid(output_orient)==0 )
  {
    printf("Error: %s is not a valid orientation code, aborting ...\n",output_orient);
    exit(0);
  }

  // if output orientation is not specified using -oo option, make it same as input image
  if(output_orient[0]=='\0')
  {
    stpcpy(output_orient, input_orient);
  }
  if(opt_v) printf("Output image orientation: %s\n",output_orient);

  // find TPIL with transforms the input image to PIL orientation
  if(opt_standard)
  {
    standard_PIL_transformation(inputimagepath, landmarksfilepath, input_orient, 0, TPIL);
  }
  else
  {
    new_PIL_transform(inputimagepath, landmarksfilepath, input_orient, TPIL, 0);
  }

  // PIL2OUT takes points from PIL space to output_orient space
  inversePILtransform(output_orient, PIL2OUT);

  //Tout transforms points from the input space to the output space
  multi(PIL2OUT, 4, 4,  TPIL, 4,  4, Tout);

  fp = fopen(outputtransformationpath,"w");
  if(fp==NULL) file_open_error(outputtransformationpath);
  printMatrix(Tout, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
  fclose(fp);

  if( reorient!='Y' && reorient!='y' && reorient!='N' && reorient!='n' )
  {
     printf("Error: The --reorient flag can only be Y or N, aborting ...\n");
  }

  // PIL2RAS takes points from PIL 2 RAS space
  // (x,y,z) in PIL -> (x,y,z) in RAS
  inversePILtransform("RAS", PIL2RAS);

  if( reorient == 'Y' || reorient == 'y')
  {
    if(opt_v) 
    {
      printf("Reorienting %s to %s orientation\n",outputimagepath, output_orient);
    }

    output_hdr = input_hdr;

    set_dim(output_dim, input_hdr);

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

    invT = inv4(Tout);
    if(!opt_nn)
    {
      output_image = resliceImage(input_image,input_dim,output_dim,invT,LIN);
    }
    else
    {
      output_image = resliceImage(input_image,input_dim,output_dim,invT,NEARN);
    }

    free(invT);

    output_hdr.pixdim[1]=odx; output_hdr.pixdim[2]=ody; output_hdr.pixdim[3]=odz;
    output_hdr.dim[1]=onx; output_hdr.dim[2]=ony; output_hdr.dim[3]=onz;
    output_hdr.magic[0]='n'; output_hdr.magic[1]='+'; output_hdr.magic[2]='1';
    sprintf(output_hdr.descrip,"Created by ART acpcdetect");
    save_nifti_image(outputimagepath, output_image, &output_hdr);

    //////////////////////////////////////////////////////////////////////////////////
    // This part of the code adjusts the SFORM and QFORM matrices of the output image
    //////////////////////////////////////////////////////////////////////////////////

    // NOTE: PIL orientation is a stepping stone here
    // (i,j,k) -> (x,y,z) in output_orient -> (x,y,z) in PIL -> (x,y,z) in RAS
    // = (i,j,k) -> (x,y,z) in RAS
    // Ttmp = PIL2RAS * OUT2PIL * Tijk2xyx

    // (i,j,k) -> (x,y,z) in output_orient
    ijk2xyz(Tijk2xyz, onx, ony, onz, odx, ody, odz);

    // OUT2PIL takes points from output_orient space to PIL space
    // (x,y,z) in output_orient -> (x,y,z) in PIL -> (x,y,z) in RAS
    PILtransform(output_orient, OUT2PIL);

    // Ttmp = PIL2RAS * OUT2PIL 
    multi(PIL2RAS, 4, 4,  OUT2PIL, 4,  4, Ttmp);

    // Ttmp = PIL2RAS * OUT2PIL * Tijk2xyx
    multi(Ttmp, 4, 4,  Tijk2xyz, 4,  4, Ttmp);

    opt_sform=YES; opt_qform=YES;
    update_qsform( (const char *)outputimagepath, Ttmp);
    //////////////////////////////////////////////////////////////////////////////////
    
    delete output_image;
  }
  else
  {
    output_hdr = input_hdr;
    output_image=input_image;

    sprintf(output_hdr.descrip,"Created by ART acpcdetect");
    save_nifti_image(outputimagepath, output_image, &output_hdr);

    // Tijk2xyz: (i,j,k) -> (x,y,z) in output_orient
    ijk2xyz(Tijk2xyz,input_dim.nx,input_dim.ny,input_dim.nz,input_dim.dx,input_dim.dy,input_dim.dz);

    // Ttmp = TPIL * Tijk2xyz : (i,j,k) -> (x,y,z) in PIL
    multi(TPIL, 4, 4,  Tijk2xyz, 4,  4, Ttmp);

    // Ttmp = PIL2RAS * TPIL * Tijk2xyz : (i,j,k) -> (x,y,z) in RAS
    multi(PIL2RAS, 4, 4,  Ttmp, 4,  4, Ttmp);

    opt_sform=YES; opt_qform=YES;
    update_qsform( (const char *)outputimagepath, Ttmp);
  }

  delete input_image;
}
