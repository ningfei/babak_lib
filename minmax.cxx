#define _minmax

int minmax(short *a, short *msk, int n, short &min, short &max)
{
   int nbv=0;

   for(int i=0; i<n; i++)
   {
      if( msk[i]>0)
      {
         min=max=a[i];
         break;
      }
   }

   for(int i=0; i<n; i++)
   {
      if( msk[i]>0)
      {
         nbv++;

         if(a[i]<min) 
         {
            min=a[i];
         }

         if(a[i]>max) 
         {
            max=a[i];
         }
      }
   }

   return(nbv);
}

