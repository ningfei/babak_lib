#ifndef _minmax_h

template<class TYPE> void minmax(TYPE *a, int n, TYPE &min, TYPE &max)
{

   min=max=a[0];

   for(int i=0; i<n; i++)
   {
      if(a[i]<min) 
      {
         min=a[i];
      }
      else if(a[i]>max) 
      {
         max=a[i];
      }
   }

   return;
}

#define _minmax_h
#endif
