#include <permutation.h>
#include <stdlib.h>		// srand48()
#include <stdio.h>
#include <time.h>		// time()

// destructor
PERMUTATION::~PERMUTATION()
{
	delete [] vec;
	delete [] ra;
}

// constructor
PERMUTATION::PERMUTATION(int s)
{
	vec = new int[s];
	ra = new double[s];

	size = s;

	for(int i=0; i<s; i++) vec[i]=i;

	srand48(time(NULL));  //	Initializes the drand function
}

void PERMUTATION::permute()
{
	for(int i=0; i<size; i++) 
	{
		vec[i]=i;
		ra[i]= drand48();        // Generates 'size' random numbers
	}

    // Heap Sort is used to sort the random numbers, resulting in permutation  of array 'vec'

	for(int i=0;i<size; i++) heap(i);

	heapsort(size);
}

void PERMUTATION::heap(int newnode)
{
  	int done=0, l;
  	double temp1;
  	int temp2;
  	l = ( newnode - 1 ) / 2;

	while( l >= 0 && !done)
  	{
    	if(ra[newnode] > ra[l])
		{
			temp1 = ra[newnode];
       		temp2 = vec[newnode];

       		ra[newnode]=ra[l];
       		vec[newnode]=vec[l];

       		ra[l]=temp1;
       		vec[l]=temp2;

       		newnode = l;
       		l = l/2;
     	}
     	else
     	done = 1;
	}
}

void PERMUTATION::heapsort(int last)
{
   	double temp1;
  	int temp2;
  
	last--;

   	while(last > 0 )
   	{
   		temp1 = ra[0];
   		temp2 = vec[0];

   		ra[0]=ra[last];
   		vec[0]=vec[last];

   		ra[last]= temp1;
   		vec[last]= temp2;

   		for(int i=0; i<last; i++)
   		heap(i);
   		last--;
   	}
}
