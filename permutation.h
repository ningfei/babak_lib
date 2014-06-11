// This is the HEADER FILE permutation.h. This is The INTERFACE for the class
// PERMUTATION, which is an ADT for 3D medical images.

#ifndef PERMUTATION_H 
#define PERMUTATION_H

class PERMUTATION
{
public:
	// constructor
	PERMUTATION(int size);	

	// destructor
	~PERMUTATION();	

	// permutation of set {0,1,2,...,size-1}
	int *vec;
	int size;

	void permute();

private:

	// array of random numbers uniformly distributed between [0.0, 1.0)
	double *ra;

	void heap(int newnode);
	void heapsort(int last);

};

#endif // PERMUTATION_H

