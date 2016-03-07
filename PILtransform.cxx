#include <babak_lib.h>
#include <sph.h>
#include <ctype.h>
#include <landmarks.h>

void transform_P(float *P, int nl, float *T)
{
   float dum[4];

   for(int i=0; i<nl; i++)
   {
      dum[0]=P[0*nl + i];
      dum[1]=P[1*nl + i];
      dum[2]=P[2*nl + i];
      dum[3]=1.0;
      multi(T,4,4,dum,4,1,dum);
      P[0*nl + i]=dum[0];
      P[1*nl + i]=dum[1];
      P[2*nl + i]=dum[2];
   }

   return;
}

// This function makes a ppm file from the MSP and displays 'nl' detected landmarks.
//im: the input image 
//lmx and lmy: are a pointer to a nl size arrays including the (x,y) for the nl landmarks
//ppmfile: the output file name
void mspPPM(SHORTIM im, int *ii, int *jj, int nl, const char *ppmfile)
{
   int yellowflag;
//   unsigned char yellow[3]={0xFF,0xFF,0x00};
   unsigned char yellow[3]={0x00,0xFF,0x00};
   unsigned char *imgTemp;
   FILE *fp;

   int kk=im.nz/2;
   
   /*size of the mark*/
   int d=4; 

   float temp;

   imgTemp=(unsigned char *)calloc(im.nv,sizeof(unsigned char));

   int low, high;
   setLowHigh(im.v, im.nv, &low, &high);

   for(int i=0;i<im.nv;i++)
   {
      if( im.v[i] >= high ) temp=255.0;
      else temp=im.v[i]*255.0/high;

      imgTemp[i]=(unsigned char)temp;
   }

   fp=fopen(ppmfile,"w");

   /*Write the header part of the PPM file*/
   fprintf(fp,"P6\n");
   fprintf(fp,"%d %d\n", im.nx, im.ny);
   fprintf(fp,"255\n");


   for (int j=0;j<im.ny;j++)
   {
      for(int i=0;i<im.nx;i++)
      {
         for(int m=0;m<nl;m++)
         {
            
            if( (i==ii[m] && jj[m]-d<j && j<jj[m]+d) || (j==jj[m] && ii[m]-d<i && i<ii[m]+d) )
            {
               fwrite(yellow,1,3,fp);
               yellowflag=1;
            }
         }
         if(yellowflag==0) 
         {
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
            fwrite(imgTemp+im.np*kk+im.nx*j+i,1,1,fp);
         }
         yellowflag=0;
      }
   }

   fclose(fp);

   delete imgTemp;
}

void convert_to_xyz(float *P, int n, SHORTIM im)
{
   /////////////////////////////////////////////////////////
   // convert P from (i,j,k) to (x,y,z) coordinates
   float T[16];
   float dum[4]={0.0,0.0,0.0,1.0};

   DIM dim;
   set_dim(dim,im);
   ijk2xyz(T, dim);

   for(int i=0; i<n; i++)
   {
      // note that P is (3 x n)
      dum[0]=P[0*n + i]; 
      dum[1]=P[1*n + i]; 
      dum[2]=P[2*n + i];
      multi(T,4,4,dum,4,1,dum);
      P[0*n + i]=dum[0]; 
      P[1*n + i]=dum[1]; 
      P[2*n + i]=dum[2];
   }
}

void convert_to_ijk(float *P, int n, SHORTIM im)
{
   /////////////////////////////////////////////////////////
   // convert P from (i,j,k) to (x,y,z) coordinates
   float T[16];
   float dum[4]={0.0,0.0,0.0,1.0};

   DIM dim;
   set_dim(dim,im);
   xyz2ijk(T, dim);

   for(int i=0; i<n; i++)
   {
      // note that P is (3 x n)
      dum[0]=P[0*n + i]; 
      dum[1]=P[1*n + i]; 
      dum[2]=P[2*n + i];
      multi(T,4,4,dum,4,1,dum);
      P[0*n + i]=dum[0]; 
      P[1*n + i]=dum[1]; 
      P[2*n + i]=dum[2];
   }
}

// matches P2=(x2,y2) as closely as possible to P1=(x1,y1) 
void point_match(float *x1, float *y1, float *x2, float *y2, int N, float *T)
{
   float params[3];
   float meanx1=0,meanx2=0;
   float meany1=0,meany2=0;
   float num=0,den=0;
   float *u1;
   float *v1;
   float *u2;
   float *v2;

   int i;

   set_to_I(T,4);

   for(i=0;i<N;i++)
   {
      meanx1+=x1[i];
      meany1+=y1[i];
      meanx2+=x2[i];
      meany2+=y2[i];
   }

   meanx1/=N;
   meany1/=N;
   meanx2/=N;
   meany2/=N;

   u1=(float *)calloc(N,sizeof(float));
   v1=(float *)calloc(N,sizeof(float));
   u2=(float *)calloc(N,sizeof(float));
   v2=(float *)calloc(N,sizeof(float));
   for(i=0;i<N;i++)
   {
      u1[i]=x1[i]-meanx1;
      v1[i]=y1[i]-meany1;
      u2[i]=x2[i]-meanx2;
      v2[i]=y2[i]-meany2;
   }

   for(i=0;i<N;i++)
   {
      num+=v1[i]*u2[i]-u1[i]*v2[i];
      den+=u1[i]*u2[i]+v1[i]*v2[i];
   }

   params[0]=atan2(num,den);
   params[1]=meanx1-cos(params[0])*meanx2+sin(params[0])*meany2;
   params[2]=meany1-sin(params[0])*meanx2-cos(params[0])*meany2;

   T[0] = cosf(params[0]);
   T[1] = -sinf(params[0]);
   T[4] = sinf(params[0]);
   T[5] = cosf(params[0]);
   T[3]=params[1];
   T[7]=params[2];

   //printMatrix(T,4,4,"T:",NULL);

   free(u1);
   free(v1);
   free(u2);
   free(v2);

   return;
}

// lmfile (landmarks file) - This is a text file with format:
// i1 j1 k1
// i2 j2 k2
// i3 j3 k3
// where 
// (i1, j1, k1) are the coordinates of the AC
// (i2, j2, k2) are the coordinates of the PC, and
// (i3, j3, k3) are the coordinates of the RP.
// subfile (subject file) - 3D T1W volume of type short in NIFTI format
// T - output 4x4 rigid-body transformation matrix that would transform
// subfile into a standardized PIL orientation
void new_PIL_transform(const char *subfile, const char *lmfile, float *T)
{
   int nl;
   float *P;  // (3 x nl) matrix 
   float *CM;
   float *invT;
   char filename[1024];
   SHORTIM subimPIL; 
   nifti_1_header PILbraincloud_hdr;
   // subfile without the directory structure and extension
   char subfile_prefix[1024]; 
   char modelfile[1024];

   if( niftiFilename(subfile_prefix, subfile)==0 ) exit(0);

   getARTHOME();

   // initial T using old PIL transformation
   standard_PIL_transformation(subfile, lmfile, 0, T);

   /////////////////////////////////////////////////////////
   // Read input volume from subfile
   SHORTIM subim; 
   nifti_1_header subim_hdr;  

   subim.v = (short *)read_nifti_image(subfile, &subim_hdr);

   if(subim.v==NULL)
   {
      printf("Error reading %s, aborting ...\n", subfile);
      exit(1);
   }

   set_dim(subim, subim_hdr);
   /////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////
   // Reslice subim to PIL space 
   sprintf(filename,"%s/PILbrain.nii",ARTHOME);
   PILbraincloud_hdr=read_NIFTI_hdr(filename);

   set_dim(subimPIL, PILbraincloud_hdr);
   invT = inv4(T);
   resliceImage(subim, subimPIL, invT, LIN);
   free(invT);

   // Detect 8 landmarks on the MSP on the subimPIL volume
   sprintf(modelfile,"%s/ali8.mdl",ARTHOME);
   sprintf(filename,"%s_PIL.nii",subfile_prefix);
   P=detect_landmarks(subimPIL, modelfile, nl);

   convert_to_xyz(P, nl, subimPIL);
 
   delete subimPIL.v;
   /////////////////////////////////////////////////////////

   {
      float avgLM[3];
      float *ST;
      float *S;
      float SST[9];
      double R[6];
      double lamda[3];
      double UT[9];
      float nrml[3]; // unit normal
      float d;

      ST = (float *)calloc(nl*3,sizeof(float));
      S =  (float *)calloc(3*nl,sizeof(float));

      transpose_matrix(P,3,nl,ST);

      avgRow(ST, nl, 3, avgLM);
      centerMatrixRow(ST, nl, 3, avgLM);

      transpose_matrix(ST,nl,3,S);

      multi(S,3,nl,ST,nl,3,SST);
      s3mat_to_vec(SST,R);

      s3eigenval(R,lamda);
      s3eigenvec(R, lamda, UT);
      
      if( UT[8] < 0.0 )
      {
         UT[6]*=-1.0; 
         UT[7]*=-1.0; 
         UT[8]*=-1.0; 
      }

      // set the unit normal equal to the eigen-vector associated
      // with the smallest eigen-value
      nrml[0]=UT[6]; 
      nrml[1]=UT[7]; 
      nrml[2]=UT[8]; 

      // compute the distance between the plane and the origin
      d=dot(nrml,avgLM,3);

      // IMPORTANT
      if(d<0.0)
      {
         nrml[0]*=-1.0; 
         nrml[1]*=-1.0; 
         nrml[2]*=-1.0; 
         d*=-1.0;
      }

      free(ST);
      free(S);
      /////////////////////////////////////////////////////////
      // project the landmarks on the plane nrml[0]*x+nrml[1]*y+nrml[2]*z=d
      float dum[4];
      float dd;
      for(int i=0; i<nl; i++)
      {
         dum[0]=P[0*nl + i]; 
         dum[1]=P[1*nl + i]; 
         dum[2]=P[2*nl + i];

         // distance between lanmark i and the place
         dd=dot(nrml,dum,3)-d;

         P[0*nl + i]=dum[0] - dd*nrml[0]; 
         P[1*nl + i]=dum[1] - dd*nrml[1]; 
         P[2*nl + i]=dum[2] - dd*nrml[2];
      }

      /////////////////////////////////////////////////////////

      {
         float theta;
         float c;
         float T2[16];
         float T1[16];
         float T21[16];
      
         set_to_I(T1,4);
         T1[3]=-d*nrml[0];
         T1[7]=-d*nrml[1];
         T1[11]=-d*nrml[2];

         if( nrml[2] > 0 )
         {
            c = nrml[2];
            if( c>1.0 ) c=1.0; // to prevent acos from getting into trouble
            if( c<0.0 ) c=0.0; // to prevent acos from getting into trouble
            theta = acosf(c);

            rotate(T2, theta, nrml[1], -nrml[0], 0.0);
         }
         else
         {
            c = -nrml[2];
            if( c>1.0 ) c=1.0; // to prevent acos from getting into trouble
            if( c<0.0 ) c=0.0; // to prevent acos from getting into trouble
            theta = acosf(c);

            rotate(T2, theta, -nrml[1], nrml[0], 0.0);
         }

         multi(T2,4,4,T1,4,4,T21);

         // update P, all landmarks will be on the x-y plane
         transform_P(P, nl, T21);

         // update T 
         multi(T21,4,4,T,4,4,T);
      }
   }

   // read CM
   {
      FILE *fp;
      int nl, r, R;
      int cm[3];

      fp=fopen(modelfile, "r");

      fread(&nl, sizeof(int), 1, fp);
      fread(&r, sizeof(int), 1, fp);
      fread(&R, sizeof(int), 1, fp);
      SPH refsph(r);

      CM = (float *)calloc(3*nl, sizeof(float));

      for(int i=0; i<nl; i++)
      {
         fread(cm, sizeof(int), 3, fp);
         fread(refsph.v, sizeof(float), refsph.n, fp);
         CM[0*nl + i]=cm[0];
         CM[1*nl + i]=cm[1];
         CM[2*nl + i]=cm[2];
      }

      fclose(fp);

      convert_to_xyz(CM, nl, subimPIL);
   }

   {
      float *x1,*x2,*y1,*y2;
      float T3[16];
      float sum;

      x1=(float *)calloc(nl,sizeof(float));
      y1=(float *)calloc(nl,sizeof(float));
      x2=(float *)calloc(nl,sizeof(float));
      y2=(float *)calloc(nl,sizeof(float));

      sum=0.0;
      for(int i=0; i<nl; i++) 
      {
         x1[i]=CM[0*nl + i];
         y1[i]=CM[1*nl + i];

         x2[i]=P[0*nl + i];
         y2[i]=P[1*nl + i];

         sum += (x1[i]-x2[i])*(x1[i]-x2[i]);
         sum += (y1[i]-y2[i])*(y1[i]-y2[i]);
      }
      //printf("%f\n",sum);

      point_match(x1, y1, x2, y2, nl, T3);

      // final T
      multi(T3,4,4,T,4,4,T);

      // update P, landmarks will be near their expected positions
      transform_P(P, nl, T3);

      //sum=0.0;
      //for(int i=0; i<nl; i++) 
      //{
      //   x1[i]=CM[0*nl + i];
      //   y1[i]=CM[1*nl + i];

      //   x2[i]=P[0*nl + i];
      //   y2[i]=P[1*nl + i];

      //   sum += (x1[i]-x2[i])*(x1[i]-x2[i]);
      //   sum += (y1[i]-y2[i])*(y1[i]-y2[i]);
      //}
      //printf("%f\n",sum);

      free(x1); free(x2); free(y1), free(y2);
   }

   // save the PIL transformation in <subfile_prefix>_PIL.mrx
   {
      FILE *fp;
      sprintf(filename,"%s_PIL.mrx",subfile_prefix);
      fp=fopen(filename,"w");
      printMatrix(T,4,4,"",fp);
      fclose(fp);
   }

   // create the *LM.ppm image
   {
      int *lmx, *lmy;

      invT = inv4(T);
      resliceImage(subim, subimPIL, invT, LIN);
      free(invT);

      //sprintf(filename,"%s_PIL2.nii",subfile_prefix);
      //save_nifti_image(filename, subimPIL.v, &PILbraincloud_hdr);

      lmx = (int *)calloc(nl,sizeof(int));
      lmy = (int *)calloc(nl,sizeof(int));
      convert_to_ijk(P, nl, subimPIL);
      for(int i=0; i<nl; i++)
      {
         lmx[i]=(int)( P[0*nl + i] + 0.5 );
         lmy[i]=(int)( P[1*nl + i] + 0.5 );
      }
      convert_to_xyz(P, nl, subimPIL);

      sprintf(filename,"%s_LM.ppm",subfile_prefix);
      mspPPM(subimPIL, lmx, lmy, nl, filename);

      delete lmx;
      delete lmy;
      delete subimPIL.v;
   }

   delete subim.v;
}

void standard_PIL_transformation(const char *imfile, const char *lmfile, int verbose, float *TPIL)
{
   char modelfile[1024]="";

   // searchradius[0] is for RP
   // searchradius[1] is for AC
   // searchradius[2] is for PC
   double searchradius[3]; // in units of mm

   // It is very important to have these initializations
   searchradius[0] = 50.0;
   searchradius[1] = 15.0;
   searchradius[2] = 15.0;

   DIM dim;
   char orient[4]="";  // empty means that the orientation is read from header
   float ac[4]={0.0, 0.0, 0.0, 1.0};
   float pc[4]={0.0, 0.0, 0.0, 1.0};
   float rp[4]={0.0, 0.0, 0.0, 1.0};

   //////////////////////////////////////////////////////////////////////////////////
   // Detect the MSP plane 
   //////////////////////////////////////////////////////////////////////////////////
   float Tmsp[16]; // takes the image to PIL orientation -- no AC/PC alignment

   findMSP(imfile,orient,lmfile,Tmsp,verbose,dim);

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // AC/PC detection 
   ///////////////////////////////////////////////////////////////////////////////////////////////

   opt_AC=opt_PC=opt_RP=YES;  // by default find these automatically

   if( lmfile[0]!='\0' )
   {
      FILE *fp;
      
      fp = fopen(lmfile,"r");
      fscanf(fp,"%f %f %f\n", &ac[0], &ac[1], &ac[2]); ac[3]=1;
      fscanf(fp,"%f %f %f\n", &pc[0], &pc[1], &pc[2]); pc[3]=1;
      fscanf(fp,"%f %f %f\n", &rp[0], &rp[1], &rp[2]); rp[3]=1;
      fclose(fp);

      opt_AC=opt_PC=opt_RP=NO;  // Do not find AC, PC, RP
   }

   opt_MSP=NO;
   detect_AC_PC_MSP(imfile,orient,modelfile,searchradius,ac,pc,rp,Tmsp,0,verbose,0);

   // convert the AC/PC from (i,j,k) in original space to (x,y,z) in PIL space
   orig_ijk_to_pil_xyz(Tmsp, dim, ac, pc);

   // 0 put the FOV center at halfway btw AC and PC to match $ARTHOME/PILbrain.nii
   ACPCtransform(TPIL, Tmsp, ac, pc, 0);

   return;
}

void PILtransform(const char *inputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<16; i++) T[i]=0.0;
   T[15]=1.0;

   for(int j=0; j<3; j++)
   {
      // read the jth element of the inputOrientCode vector
      c = inputOrientCode[j];		

      c=toupper(c);	// convert c to upper case

      // set the jth column of the T
      if( c=='P' )		T[j]   =  1.0;
      else if( c=='A' )	T[j]   = -1.0;
      else if( c=='I' )	T[4+j] =  1.0;
      else if( c=='S' )	T[4+j] = -1.0;
      else if( c=='L' )	T[8+j] =  1.0;
      else if( c=='R' )	T[8+j] = -1.0;
   }

   return;
}

void inversePILtransform(const char *outputOrientCode, float *T)
{
   char c;

   // Initialize T
   for(int i=0; i<16; i++) T[i]=0.0;
   T[15]=1.0;

   for(int i=0; i<3; i++)
   {
      // read the ith element of the outputOrientCode vector
      c = outputOrientCode[i];		
   
      c=toupper(c);	// convert c to upper case

      // set the ith row of the T
      if( c=='P' )		T[4*i] 		=  1.0;
      else if( c=='A' )	T[4*i] 		= -1.0;
      else if( c=='I' )	T[4*i + 1]	=  1.0;
      else if( c=='S' )	T[4*i + 1]	= -1.0;
      else if( c=='L' )	T[4*i + 2]	=  1.0;
      else if( c=='R' )	T[4*i + 2]	= -1.0;
   }

   return;
}
