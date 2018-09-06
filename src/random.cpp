#define _random

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()

#ifdef __MINGW32__
  #define srand48(x) srand((unsigned)(x))
  #define drand48() (rand()/(RAND_MAX + 1.0))
#endif

void initializeRandomNumberGenerator()
{
	time_t  random_seed;    // for random number generation

	// time() returns the value of time in seconds  since  00:00:00 UTC, January 1, 1970.
	// This number is used as a seed for random number generation.
	random_seed = time(NULL);

	// should be invoked before drand48()
	srand48((long)random_seed);

	// Prints today's date and time
	//printf("%s",ctime(&random_seed));
}

void sampleWithReplacementIndex(int *I, int N)
{
	for(int i=0; i<N; i++)
	{
		// generates a random interger from [0,N-1]
		I[i] = (int)(drand48()*N);
	}
}
