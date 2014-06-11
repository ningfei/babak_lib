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
#include "babak_lib.h"
#include <f2c.h>


int main(int argc, char **argv)
{
   float y[3]={0.0, 4.0, 3.0};

   removeVectorMean(y, 3);
   normalizeVector( y, 3);

   printf("y = {%f, %f, %f}\n",y[0],y[1],y[2]);
}
