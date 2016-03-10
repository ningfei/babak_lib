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
#include "../include/babak_lib.h"
#include <f2c.h>


int main(int argc, char **argv)
{
   float y[3]={1.0, 5.0, 3.0};

   double mean;

   mean = removeVectorMean( y, 3);

   printf("mean = %lf\n",mean);
   printf("y = {%f, %f, %f}\n",y[0],y[1],y[2]);
}
