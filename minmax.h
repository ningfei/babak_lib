#ifndef _minmax_h

//Input: 'array' of size 'n'
//Ouputs: The minimum and maximum values of 'array' are returned 
//in 'min' and 'max' variables.
template<class TYPE> void minmax(TYPE *array, int n, TYPE &min, TYPE &max)
{
   if(array==NULL) return;

   min=max=array[0];

   for(int i=0; i<n; i++)
   {
      if(array[i]<min) 
      {
         min=array[i];
      }
      else if(array[i]>max) 
      {
         max=array[i];
      }
   }

   return;
}

#define _minmax_h
#endif
