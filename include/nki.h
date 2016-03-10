/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  program: nki.h                                            *
*  Copyright 2004 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _nki
 
#define _nki

struct nki 
{
   char id[3];
   short flg;
   short hs;
   int nx;
   int ny;
   int nz;
   int nt;
   double dx;
   double dy;
   double dz;
   double dt;
   double cntrv[3];
   double nrmlv[3];
   double rowv[3];
   double colv[3];
   char	datatype;
   char	u;
   char	c;
   char	*cim;
   short *sim;
   int *iim;
   float *fim;
   double *dim;
};

typedef struct nki nki;

/*******************Function prototypes**************************************/
// Returns 1 if file can be read and is NKI format, 0 otherwise
int isNKI(char *file);

// returns 0 on failure, 1 on success
int saveNKI(char *filename, nki image);

// returns 0 on failure, 1 on success
int readNKI(char *filename, nki *image);
/****************************************************************************/

#endif
