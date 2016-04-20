#include <stdio.h>
#include <stdlib.h>

//*********************************************************************
// "hist1D_plot" plots the histogram of the data1 and data2 arrays based 
// on the bins. The program generates 3 files in the running directory 
// in the following formats:  name.mat 
//                            name.plt 
//                            name.png
//
// Inputs: 
//         name:  name of the output files
//         n:     number of the bins
//         bin:   array of the bin numbers
//         data1: array of the data1 values
//         data2: array of the data2 values
//
void hist1D_plot(const char *name, int n, int *bin, float *data1, float *data2)
{
   FILE *fp;
   char filename[1024]="";  // a generic filename for reading/writing stuff
   char command[1024]="";   // a generic command for running in the bash

   //////////////////////////////////////////////////////////////////////
   //Writing the data file
   {
      sprintf(filename,"%s.dat",name);
      fp = fopen(filename,"w");

      for(int i=0; i<n; i++)
         fprintf(fp,"%d %f %f\n", bin[i], data1[i], data2[i]);

      fclose(fp);
   }

   //////////////////////////////////////////////////////////////////////
   //Writing the plot file
   {
      sprintf(filename,"%s.plt",name);
      fp = fopen(filename,"w");

      fprintf(fp,"set terminal png\n");
      fprintf(fp,"set output '%s.png'\n",name);

      fprintf(fp,"set style line 1 linecolor rgb '#00ff00' lt 1 lw 1\n");
      fprintf(fp,"set style line 2 linecolor rgb '#000080' lt 3 lw 1\n");

      fprintf(fp,"plot \"%s.dat\" using 1:2 title 'Original' with lines ls 1,\\\n",name);
      fprintf(fp,"     \"%s.dat\" using 1:3 title 'Smoothed' with lines ls 2\n",name);

      fclose(fp);

   }    

   //Running the gnuplot to generate the .png file.
   sprintf(command,"gnuplot %s.plt",name);
   system(command);

}
