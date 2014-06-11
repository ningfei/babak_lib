#ifndef _sph_h

class SPH 
{
   public:
   int *i; // i voxel coordinates
   int *j; // j voxel coordinates
   int *k; // k voxel coordinates
   int n; // size in voxels
   int r; // radius in voxels
   float *v; // values at (i,j,k) - can be set by setSphere function
   SPH(int radius);
   void reset();
   void set(SHORTIM im, int ic, int jc, int kc);
   void get(float *array);
   void zeromean();
   float mean();
   float norm();
   void normalize();
   float dot(float *u);
   ~SPH();
};

#define _sph_h

#endif
