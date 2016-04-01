//*********************************************************************
// "hist2D_line_zerocrossing" plots the histogram of the baseline and follow-up intensities
// in the ROI defined by the mask image and finds the best fitting line.
// Input images should be the same size as the mask image.
// The program generates 3 files in the running directory as following:
//   hist_%s_%s.mat 
//   hist_%s_%s.plt 
//   hist_%s_%s.png
// where %s are subitute with the baseline and follow-up file names accordingly.
//
// Inputs:
//    bfile:   the baseline image file name
//    ffile:   the follow-up image file name
//    mskfile: the mask file name
// 
// Output:
//    slope:   slope of the best fitting line
// 

double hist2D_line_zerocrossing(const char *bfile, const char *ffile, const char *mskfile)
{

   FILE *fp;
   char filename[1024]=""; // a generic filename for reading/writing stuff
   char command[1024]="";  // a generic filename for reading/writing stuff

   char bprefix[1024]="";  //baseline image prefix
   char fprefix[1024]="";  //follow-up image prefix

   double slope;           //The slope of the best fitting line
   double slope_tmp;


   if(opt_v)   printf("-------------2D histogram------------\n");
   //////////////////////////////////////////////////////////////////////
   // Read mskfile from the $ARTHOME directory
   short *msk;
   DIM msk_dim;
   nifti_1_header msk_hdr;
   sprintf(filename,"%s", mskfile);
   if(opt_v) printf("Mask image:      %s\n", filename);
   msk = (short *)read_nifti_image(filename, &msk_hdr);

   if(msk==NULL)
   {
         printf("Error reading %s, aborting ...\n", filename);
         exit(1);
   }

   set_dim(msk_dim, msk_hdr);

   //////////////////////////////////////////////////////////////////////
   // Read the baseline and follow-up images
   if( niftiFilename(bprefix, bfile)==0 )  exit(0);
   if(opt_v) printf("Baseline image:  %s\n", bfile);

   if( niftiFilename(fprefix, ffile)==0 )  exit(0);
   if(opt_v) printf("Follow-up image: %s\n", ffile);

   short* bim;          // baseline image
   short* fim;          // follow-up image

   DIM bdim;            // baseline image dimensions structure
   DIM fdim;            // follow-up image dimensions structure

   nifti_1_header bhdr; // baseline image NIFTI header
   nifti_1_header fhdr; // follow-up image NIFTI header

   bim = (short *)read_nifti_image(bfile, &bhdr);
   fim = (short *)read_nifti_image(ffile, &fhdr);

   if(bim==NULL)
   {
      printf("Error reading %s, aborting ...\n", bfile);
      exit(1);
   }
   if(fim==NULL)
   {
      printf("Error reading %s, aborting ...\n", ffile);
      exit(1);
   }

   set_dim(bdim, bhdr);
   set_dim(fdim, fhdr);
  
   
   //////////////////////////////////////////////////////////////////////
   // Check the image sizes to make sure they match the mskfile
   if(bdim.nx!=msk_dim.nx || bdim.ny!=msk_dim.ny || bdim.nz!=msk_dim.nz)
   {
      printf("Size of image %s does not match the mask image, aborting ...\n", bprefix);
      exit(1);
   }
   else if(fdim.nx!=msk_dim.nx || fdim.ny!=msk_dim.ny || fdim.nz!=msk_dim.nz)
   {
      printf("Size of image %s does not match the mask image, aborting ...\n", fprefix);
      exit(1);
   }
   else
   {
      int highb=0, highf=0, highm=0;
      int ntot,n;
      int *hist;

      double *x0, *x1, *w, *w_const;
      double u[2], d;


      int high_square;
      int rangex,rangey;

      for(int i=0; i<msk_dim.nv; i++)
      if(msk[i]>0)
      {
         if(bim[i]>highb) highb=bim[i];
         if(fim[i]>highf) highf=fim[i];
         if(msk[i]>highm) highm=msk[i];
      }

      hist=(int *)calloc((highb+1)*(highf+1), sizeof(int));

      for(int i=0; i<msk_dim.nv; i++)
      if(msk[i]>0)
      {
         hist[bim[i]+(highb+1)*fim[i]]+=msk[i];
      }

      n=0;
      for(int i=0; i<(highb+1)*(highf+1); i++)
      if(hist[i]>=highm)   n++; 
      ntot=n;
   
      for(int k=0; k<=highf+highb; k++)
      {
         n=0;
         for(int i=0; i<=k; i++)
         {
            for(int j=0; j<=k; j++)
            {
               if(i<highb && j<highf && hist[i+(highb+1)*j]>=highm) 
                  n++; 
            }
         }
       
         if(n>(998*ntot)/1000)
         {
            high_square=k-1;
            break;
         }
      }
      
      rangex=rangey=high_square;

      if (high_square==0)
      {
         rangex=highb;
         rangey=highf;
      }
      else
      {
         if (high_square>highb) rangex=highb;
         if (high_square>highf) rangey=highf;
      }

      n=(rangex+1)*(rangey+1);
      x0 = (double *)calloc(n, sizeof(double));
      x1 = (double *)calloc(n, sizeof(double));
      w = (double *)calloc(n, sizeof(double));
      w_const = (double *)calloc(n, sizeof(double));

      if(opt_v)
      {
         printf("\n");
         printf("           Intensity_range  Display_range\n",highb, rangex);
         printf("Baseline:     [0,%4d]         [0,%4d]\n",highb, rangex);
         printf("Follow-up:    [0,%4d]         [0,%4d]\n",highf, rangey);
      }

      sprintf(filename,"hist_%s_%s.mat",bprefix,fprefix);
      fp = fopen(filename,"w");
      n=0;
      for(int j=0; j<=rangey; j++)
      {
         for(int i=0; i<=rangex; i++)
         {
           fprintf(fp,"%f ",(float)hist[i+(highb+1)*j]/highm);
           x0[n]=(double)i;
           x1[n]=(double)j;
           w_const[n]=(double)hist[i+(highb+1)*j]/highm;
           w[n]=1;
         
           n++;
         }
         fprintf(fp,"\n");
      }
      fclose(fp);

      //////////////////////////////////////////////////////////////////////
      //Calculating the best fitting line
      for(int iter=0; iter<3000; iter++)
      {
         line2D_norm_zerocrossing(x0, x1, w, w_const, n, u);
         slope=(-u[0]/u[1]);            //The slope of the best fitting line

         if(fabs(slope_tmp-slope)< 1.0e-6) 
         {
            if(opt_v) printf("\nNumber of iterations: %d \nSlope: %lf\n",iter, slope);
            break;
         }
         slope_tmp=slope;
      }

      free(x0);
      free(x1);
      free(w);
      free(w_const);
      free(hist);

      //////////////////////////////////////////////////////////////////////
      //Writing the plot file
      {
         sprintf(filename,"hist_%s_%s.plt",bprefix,fprefix);
         fp = fopen(filename,"w");
   
         fprintf(fp,"unset key\n");
         fprintf(fp,"plot 'hist_%s_%s.mat' matrix with image\n",bprefix,fprefix);

         fprintf(fp,"set lmargin screen 0.2\n");
         fprintf(fp,"set rmargin screen 0.9\n");
         fprintf(fp,"set tmargin screen 0.9\n");
         fprintf(fp,"set bmargin screen 0.15\n");

         fprintf(fp,"set xlabel \"Baseline intensity\"\n");
         fprintf(fp,"set ylabel \"Follow-up intensity\"\n");
         fprintf(fp,"set title \"Number of voxels\"\n");

         fprintf(fp,"if (GPVAL_DATA_Y_MAX+1 > GPVAL_DATA_X_MAX+1)\\\n");
         fprintf(fp,"    MAX=GPVAL_DATA_Y_MAX; \\\n");
         fprintf(fp,"else \\\n");
         fprintf(fp,"    MAX=GPVAL_DATA_X_MAX\n");
   
         if(rangex>rangey)
         {
            fprintf(fp,"set xrange [0:GPVAL_DATA_X_MAX]\n");
            fprintf(fp,"set yrange [0:GPVAL_DATA_X_MAX]\n");
         }
         else
         {
            fprintf(fp,"set xrange [0:GPVAL_DATA_Y_MAX]\n");
            fprintf(fp,"set yrange [0:GPVAL_DATA_Y_MAX]\n");
         }
         
         fprintf(fp,"set termoption dashed\n");
         fprintf(fp,"set style line 1 linecolor rgb '#8B008B' lt 1 lw 1\n");
         fprintf(fp,"set style line 2 linecolor rgb '#191970' lt 3 lw 1\n");
         fprintf(fp,"set style line 3 linecolor rgb '#00FF00' lt 3 lw 1\n");
         fprintf(fp,"set style line 4 linecolor rgb '#006400' lt 3 lw 1\n");

         fprintf(fp,"f(x) = (%lf - (%lf*x))/(%lf)\n",d,u[0],u[1]);
         fprintf(fp,"g(x) = x\n");
         
         fprintf(fp,"plot 'hist_%s_%s.mat' matrix with image, f(x) ls 1, g(x) ls 2\n",bprefix,fprefix);
  
         fprintf(fp,"set grid\n");
         fprintf(fp,"set style line 11 lc rgb '#252525' lt 1\n");
         fprintf(fp,"set tics nomirror out scale 0.75\n");
         fprintf(fp,"set palette defined ( 0 '#F7FBFF',\\\n");
         fprintf(fp,"                      1   '#3288BD',\\\n");
         fprintf(fp,"                      211 '#66C2A5',\\\n");
         fprintf(fp,"                      311 '#ABDDA4',\\\n");
         fprintf(fp,"                      411 '#E6F598',\\\n");
         fprintf(fp,"                      511 '#FEE08B',\\\n");
         fprintf(fp,"                      611 '#FDAE61',\\\n");
         fprintf(fp,"                      711 '#F46D43',\\\n");
         fprintf(fp,"                      811 '#D53E4F')\n");

         fprintf(fp,"set size ratio -1\n");
         fprintf(fp,"set terminal png font 'arial,28' enhanced nocrop size 2000,1600\n");
         
         fprintf(fp,"set output 'hist_%s_%s.png'\n",bprefix,fprefix);
         fprintf(fp,"replot\n");
         fclose(fp);
      }    

      //Checking the .plt file
      sprintf(filename,"hist_%s_%s.plt",bprefix,fprefix);
      if(!checkFileExistence(filename) && opt_v)
         printf("Error in writing the %s file.\n",filename);

      //Running the gnuplot to generate the .png file.
      sprintf(command,"gnuplot hist_%s_%s.plt",bprefix,fprefix);
      system(command);

      //Checking the .png file
      sprintf(filename,"hist_%s_%s.png",bprefix,fprefix);
      if(!checkFileExistence(filename) && opt_v)
         printf("Error in running the gnuplot.\n File %s is not generated.",filename);
   }

   if(opt_v)   printf("-----------------------------------\n");
   return slope;
}
