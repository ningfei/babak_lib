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
   {"--do-not-reorient", 0, 'R'},
   {"-lm", 1, 'L'},
   {"-landmarks", 1, 'L'},
   {"--landmarks", 1, 'L'},
   {"--center-AC", 0, 'M'},
   {"-center-AC", 0, 'M'},
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
   {"-input-orient",1,'O'}, // secret option
   {"--input-orient",1,'O'}, // secret option
   {"-io",1,'O'}, // secret option
   {"-h",0,'h'},
   {"--help",0,'h'},
   {"-help",0,'h'},
   {"-o",1,'o'},
   {"-output",1,'o'},
   {"--output",1,'o'},
   {"-oo",1,'u'},
   {"-output-orient",1,'u'},
   {"--output-orient",1,'u'},
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
   "[-o/-output <output volume> -oo <output orientation code>]\n"
   "[-nx <int> -ny <int> -nz <int> -dx <float> -dy <float> -dz <float> -noppm -notxt]\n"
   "[--do-not-reorient]\n"
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
   "-o or -output <filename>: If this option is present, the program outputs an AC/PC and MSP aligned image\n"
   "with the given filename.\n"
   "-oo <output orientation code>: Three-letter orientation code of the output image.  If this is not\n"
   "specified, then the output image will have the same orientation as the input image. Ex:\n"
   "\t\tPIL for Posterior-Inferior-Left\n"
   "\t\tRAS for Right-Anterior-Superior\n"
   "-nx, -ny, -nz: x/y/z matrix dimensions of the output image (default=same as input image).\n"
   "-dx, -dy, -dz: x/y/z voxel dimensions of the output image (default=same as input image).\n"
   "--center-AC: Make AC the center of the output FOV.\n\n"
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
  char opt_reorient=YES;
  int nPA, nLR, nIS;
  float dPA, dLR, dIS;
  char opt_standard=NO;
  FILE *fp;
  short *ipimage;
  DIM ipdim, opdim;  // input (ip) and output (op) dimension structures (dim)
  nifti_1_header iphdr, ophdr; // 348 bytes
  float PIL2RAS[16];
//  int opt_T2=NO;
  char optransformationpath[1024]="";
  char landmarksfilepath[512]="";
//  char modelfile[1024]="";
  char ipimagepath[1024]="";
  char ipimagename[512]="";
  char ipimagedir[512]="";  // important to initialize to ""
  char opimagepath[1024]="";
  char iporient[4]="";
  char oporient[4]="";
  float Tout[16]; // transforms ipimage to the specified output orientation
  float TPIL[16]; // transforms ipimage to PIL orientation

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
            opt_reorient=NO;
            break;
         case 'L':
            sprintf(landmarksfilepath,"%s",optarg);
            break;
         case 'o':
            sprintf(opimagepath,"%s",optarg);
            break;
         case 'i':
            sprintf(ipimagepath,"%s",optarg);
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
            sprintf(iporient,"%s",optarg);
            break;
         case 'u':
            sprintf(oporient,"%s",optarg);
            break;
//         case 'm':
//            sprintf(modelfile,"%s",optarg);
//            break;
         case '?':
            print_help_and_exit();
      }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // This block of code set the following variables: ipimagepath, ipimagename, ipimagedir
  // iporient, ipimage, iphdr, ipdim, nPA, nIS, nLR, dPA, dIS, dLR 

  if( ipimagepath[0]=='\0' )
  {
    printf("Please specify an input image using: -i <inputimage.nii>\n");
    exit(0);
  }
  if(opt_v) printf("Input image: %s\n",ipimagepath);

  // determine input image filename without the .nii suffix
  if( niftiFilename(ipimagename, ipimagepath)==0 ) { exit(1); }

  // determine input image directory
  getDirectoryName(ipimagepath, ipimagedir);

  // if input orientation is specified using --input-orient option, make sure it's valid
  // --input-orient overrides the input orientation inferred from the image header
  if(iporient[0]!='\0' && isOrientationCodeValid(iporient)==0)
  {
    printf("%s is not a valid orientation code, aborting ...\n",iporient);
    exit(0);
  }

  // if input orientation is not specified using --input-orient option, read it from input image
  // this is almost always going to be the case
  if(iporient[0]=='\0')
  {
    getNiftiImageOrientation(ipimagepath, iporient);
  }

  // change iporient to upper case letters, just in case
  iporient[0]=(char)toupper((int)iporient[0]);
  iporient[1]=(char)toupper((int)iporient[1]);
  iporient[2]=(char)toupper((int)iporient[2]);

  if(opt_v) printf("Input image orientation: %s\n",iporient);

  ipimage = (short  *)read_nifti_image(ipimagepath, &iphdr);
  if(ipimage==NULL)
  {
    printf("Error reading %s, aborting ...\n", ipimagepath);
    exit(1);
  }

  if( iphdr.datatype != DT_SIGNED_SHORT && iphdr.datatype != DT_UINT16)
  {
    printf("This program can only handle short and unsigned short data types, aborting ...\n");
    exit(0);
  }

  set_dim(ipdim, iphdr);

  if(opt_v) 
  {
    printf("Input image matrix size: %d x %d x %d\n",ipdim.nx,ipdim.ny,ipdim.nz);
    printf("Input image voxel size: %6.4f x %6.4f x %6.4f\n",ipdim.dx,ipdim.dy,ipdim.dz);
  }

  // set nPA, nIS, nLR, dPA, dIS, dLR depending on iporient
  if(iporient[0]=='P' || iporient[0]=='A') { nPA=ipdim.nx; dPA=ipdim.dx; }
  if(iporient[0]=='I' || iporient[0]=='S') { nIS=ipdim.nx; dIS=ipdim.dx; }
  if(iporient[0]=='L' || iporient[0]=='R') { nLR=ipdim.nx; dLR=ipdim.dx; }
  if(iporient[1]=='P' || iporient[1]=='A') { nPA=ipdim.ny; dPA=ipdim.dy; }
  if(iporient[1]=='I' || iporient[1]=='S') { nIS=ipdim.ny; dIS=ipdim.dy; }
  if(iporient[1]=='L' || iporient[1]=='R') { nLR=ipdim.ny; dLR=ipdim.dy; }
  if(iporient[2]=='P' || iporient[2]=='A') { nPA=ipdim.nz; dPA=ipdim.dz; }
  if(iporient[2]=='I' || iporient[2]=='S') { nIS=ipdim.nz; dIS=ipdim.dz; }
  if(iporient[2]=='L' || iporient[2]=='R') { nLR=ipdim.nz; dLR=ipdim.dz; }
  ///////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  // This block of code sets the following variables: oporient, opdim, opimagepath

  // if ouput orientation is specified using -oo option, make sure it's valid
  if(oporient[0]!='\0' && isOrientationCodeValid(oporient)==0 )
  {
    printf("%s is not a valid orientation code, aborting ...\n",oporient);
    exit(0);
  }

  // if output orientation is not specified using -oo option, make it same as input image
  if(oporient[0]=='\0')
  {
    stpcpy(oporient, iporient);
  }

  oporient[0]=(char)toupper((int)oporient[0]);
  oporient[1]=(char)toupper((int)oporient[1]);
  oporient[2]=(char)toupper((int)oporient[2]);

  if(oporient[0]=='P' || oporient[0]=='A') { opdim.nx=nPA; opdim.dx=dPA; }
  if(oporient[0]=='I' || oporient[0]=='S') { opdim.nx=nIS; opdim.dx=dIS; }
  if(oporient[0]=='L' || oporient[0]=='R') { opdim.nx=nLR; opdim.dx=dLR; }
  if(oporient[1]=='P' || oporient[1]=='A') { opdim.ny=nPA; opdim.dy=dPA; }
  if(oporient[1]=='I' || oporient[1]=='S') { opdim.ny=nIS; opdim.dy=dIS; }
  if(oporient[1]=='L' || oporient[1]=='R') { opdim.ny=nLR; opdim.dy=dLR; }
  if(oporient[2]=='P' || oporient[2]=='A') { opdim.nz=nPA; opdim.dz=dPA; }
  if(oporient[2]=='I' || oporient[2]=='S') { opdim.nz=nIS; opdim.dz=dIS; }
  if(oporient[2]=='L' || oporient[2]=='R') { opdim.nz=nLR; opdim.dz=dLR; }
  if(onx > 0) opdim.nx=onx; 
  if(ony > 0) opdim.ny=ony; 
  if(onz > 0) opdim.nz=onz;
  if(odx > 0.0) opdim.dx=odx; 
  if(ody > 0.0) opdim.dy=ody; 
  if(odz > 0.0) opdim.dz=odz;
  opdim.nt=1; 
  opdim.dt=0.0; 
  opdim.np=opdim.nx*opdim.ny; 
  opdim.nv=opdim.np*opdim.nz; 

  // if this option is selected, the output image size and orientation will be the same as
  // the input image, over-ridding all other settings
  if(opt_reorient==NO)
  {
    set_dim(opdim, iphdr);
    stpcpy(oporient, iporient);
  }

  if( opimagepath[0]=='\0' )
  {
    sprintf(opimagepath,"%s/%s_%s.nii",ipimagedir,ipimagename,oporient);
  }

  if(opt_v) 
  {
    printf("Output image: %s\n",opimagepath);
    printf("Output image matrix size: %d x %d x %d\n",opdim.nx,opdim.ny,opdim.nz);
    printf("Output image voxel size: %6.4f x %6.4f x %6.4f\n",opdim.dx,opdim.dy,opdim.dz);
    printf("Output image orientation: %s\n",oporient);
  }
  ///////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  if(searchradius[0]<=0 || searchradius[0]>200.0) searchradius[0]=50.0;
  if(searchradius[1]<=0 || searchradius[1]>100.0) searchradius[1]=15.0;
  if(searchradius[2]<=0 || searchradius[2]>100.0) searchradius[2]=15.0;

  if(opt_v && landmarksfilepath[0]!='\0') 
  {
    printf("Manually specified landmarks: %s\n",landmarksfilepath);
  }

  // find TPIL which transforms the input image to PIL orientation
  if(opt_standard)
  {
    standard_PIL_transformation(ipimagepath, landmarksfilepath, iporient, 0, TPIL);
  }
  else
  {
    new_PIL_transform(ipimagepath, landmarksfilepath, iporient, TPIL, 0);
  }
  ///////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  {
    float Tout_FSL[16];
    float PIL2OUT[16];

    sprintf(optransformationpath,"%s/%s.mrx",ipimagedir,ipimagename);
    if(opt_v) printf("Output transformation matrix: %s\n",optransformationpath);

    // PIL2OUT takes points from PIL space to oporient space
    inversePILtransform(oporient, PIL2OUT);

    //Tout transforms points from the input space to the output space
    multi(PIL2OUT, 4, 4,  TPIL, 4,  4, Tout);

    fp = fopen(optransformationpath,"w");
    if(fp==NULL) file_open_error(optransformationpath);
    printMatrix(Tout, 4, 4, "ART acpcdetect tilt correction matrix:", fp);
    fclose(fp);

    art_to_fsl(Tout, Tout_FSL, ipdim, opdim);
    sprintf(optransformationpath,"%s/%s_FSL.mat",ipimagedir,ipimagename);
    if(opt_v) printf("Output transformation matrix (FSL format): %s\n",optransformationpath);
    fp=fopen(optransformationpath,"w");
    if(fp==NULL) file_open_error(optransformationpath);
    printMatrix(Tout_FSL,4,4,"",fp);
    fclose(fp);
  }
  ///////////////////////////////////////////////////////////////////////////////////////

  // PIL2RAS takes points from PIL 2 RAS space
  // (x,y,z) in PIL -> (x,y,z) in RAS
  inversePILtransform("RAS", PIL2RAS);

  if( opt_reorient == YES)
  {
    short *opimage;
    float *invT;
    float Tijk2xyz[16];
    float Ttmp[16];
    float OUT2PIL[16];

    ////////////////////////////////////////////////////////////////////
    // setup ophdr
    ophdr = iphdr;
    ophdr.pixdim[1]=opdim.dx; 
    ophdr.pixdim[2]=opdim.dy; 
    ophdr.pixdim[3]=opdim.dz;
    ophdr.dim[1]=opdim.nx; 
    ophdr.dim[2]=opdim.ny; 
    ophdr.dim[3]=opdim.nz;
    sprintf(ophdr.descrip,"Created by ART acpcdetect");

    // NOTE: PIL orientation is a stepping stone here
    // (i,j,k) -> (x,y,z) in oporient -> (x,y,z) in PIL -> (x,y,z) in RAS
    // = (i,j,k) -> (x,y,z) in RAS
    // set Ttmp = PIL2RAS * OUT2PIL * Tijk2xyx

    // (i,j,k) -> (x,y,z) in oporient
    ijk2xyz(Tijk2xyz, opdim);

    // OUT2PIL takes points from oporient space to PIL space
    // (x,y,z) in oporient -> (x,y,z) in PIL -> (x,y,z) in RAS
    PILtransform(oporient, OUT2PIL);

    // set Ttmp = PIL2RAS * OUT2PIL 
    multi(PIL2RAS, 4, 4,  OUT2PIL, 4,  4, Ttmp);

    // Ttmp = PIL2RAS * OUT2PIL * Tijk2xyx
    multi(Ttmp, 4, 4,  Tijk2xyz, 4,  4, Ttmp);

    update_qsform(ophdr, Ttmp);
    ////////////////////////////////////////////////////////////////////

    invT = inv4(Tout);
    if(!opt_nn)
    {
      opimage = resliceImage(ipimage,ipdim,opdim,invT,LIN);
    }
    else
    {
      opimage = resliceImage(ipimage,ipdim,opdim,invT,NEARN);
    }
    free(invT);

    save_nifti_image(opimagepath, opimage, &ophdr);
    
    delete opimage;
  }

  if( opt_reorient == NO)
  {
    float Tijk2xyz[16];
    float Ttmp[16];

    ////////////////////////////////////////////////////////////////////
    // setup ophdr
    ophdr = iphdr;
    sprintf(ophdr.descrip,"Created by ART acpcdetect");

    // Tijk2xyz: (i,j,k) -> (x,y,z) in oporient
    ijk2xyz(Tijk2xyz, ipdim);

    // Ttmp = TPIL * Tijk2xyz : (i,j,k) -> (x,y,z) in PIL
    multi(TPIL, 4, 4,  Tijk2xyz, 4,  4, Ttmp);

    // Ttmp = PIL2RAS * TPIL * Tijk2xyz : (i,j,k) -> (x,y,z) in RAS
    multi(PIL2RAS, 4, 4,  Ttmp, 4,  4, Ttmp);

    update_qsform(ophdr, Ttmp);
    ////////////////////////////////////////////////////////////////////

    save_nifti_image(opimagepath, ipimage, &ophdr);
  }

  delete ipimage;
}
