#ifndef _babak_lib_h
#define _babak_lib_h

#include "minmax.h"
#include "stats.h"
#include "nki.h"
#include "permutation.h"
#include "volume.h"
#include <nifti1.h>
#include <nifti1_io.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef unsigned short uint2;
typedef unsigned int uint4;
typedef short int2;
typedef int int4;
typedef float float4;
typedef double float8;

//////////////////////////////////////////////////////////////////////////////////////////////////

#define YES 1
#define NO 0
#define ESMALL 1e-10

//////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct shortim
{
   int nx;
   int ny;
   int nz;
   int np;
   int nv;
   float dx;
   float dy;
   float dz;
   short *v; // image values
} SHORTIM;

struct DICOM_file_meta_info
{
   // maximum bytes for VR=UI is 64,  1 byte is added for the string terminator '\0'
   char media_storage_SOP_class[65];  // Tag (0002,0002)
   char transfer_syntax[65]; // Tag (0002,0010)
   int dataset_offset;
};
typedef struct DICOM_file_meta_info DICOM_file_meta_info;

struct DICOM_hdr
{
   char patientID[65];    // Tag (0010,0020)  VR=LO (Long String 64 chars max.)
   float slice_thickness; // Tag (0018,0050)
   int TR;             // Tag (018,0080)
   int TE;             // Tag (018,0081)
   int TI;             // Tag (018,0082)
   int NEX;            // Tag (018,0083)
   float dz;           // Tag (018,0088)
   float flip_angle;   // Tag (018,1314)
   char seriesID[65];  // Tag (020,000E)
   int series_number;  // Tag (020,0011)
   float TLHC[3];      // Tag (020,0032)
   float rowvec[3];    // Tag (020,0037)
   float colvec[3];    // Tag (020,0037)
   char frame_of_referenceID[65]; // Tag (020,0052)
   int nz; // Tag (020,1002)
   float slice_location; // Tag (020,1041)
   int ny; // Tag (028,0010)
   int nx; // Tag (028,0011)
   float dy; // Tag (028,0030)
   float dx; // Tag (028,0030)
};
typedef struct DICOM_hdr DICOM_hdr;

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
struct model_file_hdr
{
   int nxHR;
   int nzHR;
   float dxHR;
   int nxLR;
   int nvol; // number of image volumes in the training set 
   int RPtemplateradius;
   int RPtemplateheight;
   int RPtemplatesize;
   int ACtemplateradius;
   int ACtemplateheight;
   int ACtemplatesize;
   int PCtemplateradius;
   int PCtemplateheight;
   int PCtemplatesize;
   int nangles; // number of angles, each template is rotated by this many angles and saved
};

typedef struct model_file_hdr model_file_hdr;

struct model_file_tail
{
   float RPPCmean[2]; // RPPC is a vector on the MSP that points from the RP point to the PC.   RP------->PC
   float parcomMean;
   float percomMean; 
   float RPmean[2];
};

typedef struct model_file_tail model_file_tail;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
struct im_params {

   int MRimageId;
   int PETimageId;

	/* maximum values of the image 1 */	
	unsigned char max1; 	

	int 	*Size;   /* size of non-zero CC's in CCI */
	float *S,*SS;  

	int 	FROM;
	int 	TO;

	/* number of image 1 columns, rows and slices  */
	int	nx1,ny1,nz1; 

	/* number of image 1 voxels, number of image 1 pixels per slice */
	int	nv1,np1;      	

	float  dx1,dy1,dz1; 	/* image 1 voxel size (mm) */

	/* coordinates of image 1 volume center wrt pixel (0,0,0) in (mm) */
	float  xc1,yc1,zc1; 	

	// number of image 2 columns, rows and slices  
	int	nx2,ny2,nz2; 
	float  dx2,dy2,dz2; 	// image 2 voxel size (mm) 

	/* number of image 2 voxels, number of image 2 pixels per slice */
	int	nv2,np2;      

	/* coordinates of image 2 volume center wrt pixel (0,0,0) in (mm) */
	float  xc2,yc2,zc2; 	

	int NCC;

	int NP;  // number of non-zero pixels in CCI

	int sf;  // speed factor 

	unsigned *CCI;  	/* labeled connected compontents' image */

	unsigned char *data1;  /* image 1 */
	unsigned char *data2;  /* image 2 */
};
//////////////////////////////////////////////////////////////////////////////////////////////////

struct dicominfo
{
   char patientName[321];
   char DOB[11];
   char patientID[65];
   char studyDate[11];
   char TE[17];
   char TR[17];
   char TI[17];
   char ETL[13];
   char NEX[17];
   char flipAngle[17];
   char bandwidth[17];
   char freq[17];
   unsigned short acquisitionMatrix[4];
   char phaseFOV[17];
   char dum[3];
};
typedef struct dicominfo dicominfo;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct dim {
   int nx; // number of columns 
   int ny; // number of rows
   int nz; // number of slices
   int nt; // number of frames (epochs)
   float dx; // x voxel dimension (mm)
   float dy; // y voxel dimension (mm)
   float dz; // z voxel dimension (mm)
   float dt; // t between frames (sec)
} DIM;
//////////////////////////////////////////////////////////////////////////////////////////////////

// Input: (x,y,z) a vector defined in RAS system
// Output: One of six charaters {R,L,A,P,S,I}
char directionCode(float x, float y, float z);
void getNiftiImageOrientation(const char *filename, char *orientation);
int checkNiftiFileExtension(const char *filename);
void memory_allocation_error(const char *s);
void file_open_error(const char *s);
void errorMessage(const char *s);
int isOrientationCodeValid(const char *orientCode);
void swap_model_file_hdr(model_file_hdr *hdr);
void swap_model_file_tail(model_file_tail *tail);
void PILtransform(const char *orientCode, float *orientMat);
void inversePILtransform(const char *orientCode, float *orientMat);
short *reorientVolume(short *v1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1, float *orientMat,
int &nx2, int &ny2, int &nz2, float &dx2, float &dy2, float &dz2);
void rotate(float *T, float alpha, float x, float y, float z);
float *rotate(float alpha, float x, float y, float z);
void rotate(float *T, float CosAlpha, float SinAlpha, float x, float y, float z);
void setLowHigh(short *image, int nv, int *low, int *high);	
void setLowHigh(short *image, int nv, int *low, int *high, float percent);
void compute_cm(short *image, int nx, int ny, int nz, float dx, float dy, float dz, float *x, float *y, float *z);
void standardize(float *x, int n);

#ifndef _singular_value_decomposition
extern int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U, 
                      double* singular_values, double* V, double* dummy_array);
extern void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,  
                        double tolerance, int nrows, int ncols, double *Astar);
#endif

#ifndef _DKI
extern void formA7(float *A, int n, float *u1, float *u2, float *u3, float *b, float *w);

extern void restricted_tensor_model(float *x, float *y);
extern void second_and_fourth_order_moments(float *x, float *y, float *E, float *F, int v);
extern float mardia_kurtosis(float *x1, float *x2, int v);
extern float mardia_kurtosis(float *x);

extern float *formA(float *A, int n, float *u1, float *u2, float *u3, float *b, float scale);
extern float estimate_K(float *A, float *x, float *y, int n);
extern void update_x(float *A, float *x, float *y, float K, int n);
extern void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v, float scale);
extern void update_v(float *u0, float *u1, float *u2, float *b, float *y, int n, float *v1, float *v2, float *v3, float scale);
extern float quadratic_form_1(float *x, float *v);
extern float quadratic_form_1(float x0, float x1, float x2, float *v);
extern void vector_to_symmetric_matrix_form(float *u, float *D);
extern void symmetric_matrix_to_vector_form(float *D, float *u);
extern void cholesky_composition_1(float *v, float *u);
extern void cholesky_composition_2(float *L, float *D);
extern int cholesky_decomposition_1(float *u, float *v);
extern int cholesky_decomposition_2(float *D, float *L);
extern void positive_definite_tensor_model(float *v1, float *v2, float *y);
#endif

#ifndef _artlib
extern char opt_ppm;
extern char opt_txt;
extern char opt_AC;
extern char opt_PC;
extern char opt_RP;
extern void sub2trg_rigid_body_transformation(float *sub2trg, const char *subfile, const char *trgfile);
extern void getDirectoryName(const char *pathname, char *dirname);
extern void forwardTCSAP(float *xvec, float *yvec, float *zvec, float *TLHC, float *angle, float *translation, DIM dim);
extern void backwardTCSAP(float *xvec, float *yvec, float *angle);
extern void orig_ijk_to_pil_xyz(float *Tmsp, DIM orig_dim, float *AC, float *PC);
extern void initialAC(float Ax, float Ay, float Bx, float By, float *Cx, float *Cy, float parcomMean, float percomMean);
extern short *thresholdImageOtsu(short *im, int nv, int *nbv);
extern void defineTemplate(int r, int h, short *x, short *y, short *z);
extern char *defineACregion(DIM dim, float *RP, float *PC, float parcomMean, float percomMean, double ACsr);
extern char *definePCregion(DIM HR, float *RP, float *RPPCmean, double PCsr);
extern char *expandMask(short *mask_HR, DIM HR, float *RPmean, double RPsr);
extern void ACPCtransform(float *Tacpc, float *Tmsp, float *AC, float *PC, char flg);
extern void compute_MSP_parameters_from_Tmsp(float *Tmsp, float *n, float *d);
extern int detect_AC_PC_MSP( const char *imagefilename, char *orientation, char *modelfile, double *searchradius,
float *AC, float *PC, float *RP, float *Tmsp, int opt_D, int opt_v, int opt_T2);
extern float optimizeNormalVector(short *image,DIM dim,float *A, float *B, float *C);
extern float reflection_cross_correlation2(short *image, DIM dim, float A, float B, float C);
extern float reflection_cross_correlation(short *image, DIM dim, float a, float b, float c, float d);
extern void findInitialNormalVector(short *image, DIM dim, float *A, float *B,float *C);
extern float msp(short *im_in, int nx, int ny, int nz, float dx, float dy, float dz, float *A, float *B, float *C);
extern void computeTmsp(char *orientation, short *volOrig, DIM dim, float *Tmsp);
extern int save_as_ppm(const char *filename, int nx, int ny, char *R, char *G, char *B);
extern int save_as_ppm(const char *filename, int nx, int ny, unsigned char *R, unsigned char *G, unsigned char *B);
extern void combine_warps_and_trans(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
#endif

#ifndef _binomial
extern double BinomialCoeff( int n, int k );
extern void k_subset ( int k, double rank, int a[] );
#endif

#ifndef _maskOps
extern char *find_foreground_mask(short *im, int nv, int nb, int nclass, int niter, short *thresh);
extern short *mask_image(short *im, char *mask, int nv, int *nbv_out);
extern void cc_filter(char *mask, int nx, int ny, int nz);
#endif

#ifndef _EMFIT
extern void EMFIT1d(double *hist, double *fit, short *label, int nb, double *mean, double *var, double *p, int nclass, int niter);
#endif

#ifndef _max_cc

extern void Connected_Components( char *im, int nx, int ny, int nz, int *LABEL, int *N, int **Clabel, int **Size);
extern void max_Connected_Component(char *im, int nx, int ny, int nz,int *ncc, int *maxsize);
extern void thr_Connected_Component(char *im,int thr, int nx, int ny, int nz,int *ncc, int *ntcc);
extern void heightthr_Connected_Component(float *im, float thr, int nx, int ny, int nz, int *ncc, int *ntcc);
extern void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L);
extern void Connected_Component_location(char *im, int nx, int ny, int nz, int *ncc, int **S, int **L, int **CCsize);

#endif

#ifndef _statistics
extern float imageMean(short *im, short *msk, int nv);
extern double one_sample_t(double *x, int n);

extern float median(float *x, char *mask, int n);
extern void scaleAbsToOne(double *y, int n, int p);
extern double scaleAbsToOne(double *y, int n);
extern void decomposeVector(double *x, double *xpar, double *xper, double *u, int n);

extern void removeComponent(double *x, double *y, int n);
extern void removeComponent(float *x, float *y, int n);

extern double removeVectorMean(float *y, int n);
extern double removeVectorMean(double *y, int n);
extern double removeVectorMean(float *x, float *y, int n);
extern double removeVectorMean(double *x, double *y, int n);
extern double removeVectorMean(short *x, double *y, int n);
extern void removeVectorMean(double *y, int n, int p);

extern void partialCorrelation(double *Y, double *X1, double *X2, int n, double *pr1, double *pr2);
extern void partialCorrelation(float *Y, float *X1, float *X2, int n, double *pr1, double *pr2);
extern void partialCorrelation(short *Y, float *X1, float *X2, int n, double *pr1, double *pr2);

#endif

#ifndef _ginverse
extern int ginverse(float *X, int N, int p, float *G);
extern int ginverse(double *X, int N, int p, float *G);
#endif

#ifndef _convolution
extern float conv_pnt_sk(unsigned char *x,int sx,float *h,int sh,int i0);
extern float conv_pnt_sk(short *x,int sx,float *h,int sh,int i0);
extern float conv_pnt_sk(float *x,int sx,float *h,int sh,int i0);
extern float *conv_sk(short *x,int sx,float *h,int sh);
extern float *conv_sk(float *x,int sx,float *h,int sh);
extern void conv_sk(float *x, float *y, int sx,float *h,int sh);
extern void conv_sk_inplace(float *x, int sx, float *h, int sh);
#endif

#ifndef _gaussian_kernel
extern float *gaussian_kernel(float sd, int *n);
#endif

#ifndef _hpsort
extern void hpsort(unsigned long n, float *ra);
extern void hpsort(unsigned long n, float *ra, int *indx);
#endif

#ifndef _random
extern void initializeRandomNumberGenerator();
extern void sampleWithReplacementIndex(int *I, int N);
#endif

#ifndef _dicomIO

#define IMPLICIT_LITTLE_ENDIAN 0
#define EXPLICIT_LITTLE_ENDIAN 1
#define EXPLICIT_BIG_ENDIAN 2

extern int readFileMetaInfo(const char *filename, DICOM_file_meta_info *file_meta_info, char v);
extern int readDataSet(const char *filename, DICOM_hdr *hdr, char v);

extern void readDicomInfo(const char *file, int np, dicominfo *info);
extern int readImageParams(const char *file, float *TLHC, float *rowvec, float *colvec, 
float *dx, float *dy, float *dz, int *TE, char *patientID, int *imageNumber, int *seriesNumber, int np);
extern int readPhaseEncodingDirection(const char *file, char *PED, int np);
extern int readPixelData(const char *file, char *data, int opt_j, int np);
extern int readTLHC(const char *file, float *TLHC, int np);
extern int readRowCol(const char *file, float *rowvec, float *colvec, int np);
extern int readNEX(const char *file, int *nex);
extern int readSliceThickness(const char *file, float *sliceThickness);
extern int readSeriesNumber(const char *file, int *seriesNumber, int np);
extern int readImageNumber(const char *file, int *imageNumber, int np);
extern int readTE(const char *file, int *TE, int np);
extern int readTR(const char *file, int *TR);
extern int readImageVoxelSize(const char *file, float *dx, float *dy, float *dz, int np);
extern int readImageMatrixSize(const char *file, unsigned short *nx, unsigned short *ny);
extern int readPatientID(const char *file, char *patientID, int np);
extern int dicomFormat(const char *file);
extern int readMetaFileInfo(const char *file, long *offset, int *syntax);
extern int readTag(const char *file, unsigned short iGN, unsigned short iEN, long byteOffset, int transferSyntax,
unsigned long *oVL, char *oV, long *valueOffset);
extern int read_element(char *filename, short S_GN, short S_EN, char *V, int *VL);
extern void readMatrixSize(char *filename, int *nx, int *ny);
extern void readVoxelSize(char *filename, float *dx, float *dy, float *dz);
extern int readImageSliceThickness(const char *file, float *dz, int np);

#endif

#ifndef _nkiIO
extern int isNKI(char *file);
extern int saveNKI(char *filename, nki image);
extern int readNKI(char *filename, nki *image);
#endif

#ifndef _subsets
extern int *subsets(int N, int M);
extern int binomialCoeff(int N,int M);
#endif

#ifndef _cubicspline
extern float cubicSplineSynthesis(float *c, int nx, int ny, int nz, float x, float y, float z, float *beta, float del);

extern short *resliceImageCubicSpline(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *resliceImageCubicSpline(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *computeBeta(float *del);

extern void cubicSplineAnalysis(unsigned char *s, float *c, int nx, int ny, int nz);
extern void cubicSplineAnalysis(short *s, float *c, int nx, int ny, int nz);
extern void cubicSplineAnalysis(float *s, float *c, int nx, int ny, int nz);

extern void cubicSplineAnalysis(float *s, float *c, int N);

#endif

#ifndef _registration
extern void set_transformation(float x, float y, float z, float ax, float ay, float az, const char *code, float *T);
extern int label_CCI(short *KMI, int size_thresh,struct im_params * IP, int nvoxels);
extern short (*interpolator)(float x, float y, float z, short *array, int nx, int ny, int nz, int np);
extern float P[12];
extern struct im_params IP;

extern unsigned char linearInterpolatorUC(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np);
extern void scale_short_minmax(short *imagein, unsigned char **imageout, int np, int min,int max);

extern void testCostFunc1(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz,
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz);

extern float *transformation(float x, float y, float z, float ax, float ay,
float az, float sx, float sy, float sz, int rX, int rY, int rZ, char *code);

extern void cca(short *im, int nx, int ny, int nz);

extern int loadTransformation( char *filename, float *T);
extern int findThresholdLevel(short *image_in, int nv);
extern float *findTransMatrix(short *trg, int Tnx, int Tny, int Tnz, float Tdx, float Tdy, float Tdz,
short *obj, int Onx, int Ony, int Onz, float Odx, float Ody, float Odz);

extern short *KMcluster(short *ccImage, short *im_in, int nclass, int maxiter, int thresh, int Low, int High, int nv);
#endif

#ifndef _legendre
//extern float *FourierLegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz);
//extern double *FourierLegendreAnalysis(float *f, int nx, int ny, int nz, int mx, int my, int mz, int N);
extern double *LegendreAnalysis(float *image, int nx, int ny, int nz, int mx, int my, int mz);
extern	double *LegendreAnalysis(float *image, int nx, int ny, int mx, int my);
extern void LegendreSynthesis(double *c, int mx, int my, float *image, int nx, int ny);
extern float *LegendreSynthesis(double *c, int nx, int ny, int nz, int mx, int my, int mz);
void LegendrePoly(double *p0, double *q0, double *p1, double *q1, double *x, int N, int n);
void LegendrePoly(float *p0, float *q0, float *p1, float *q1, float *x, int N, int n);
void integral_1d(double *a, short *b, int n, double *d);
void integral_1d(double *a, float *b, int n, double *d);
#endif

#ifndef _matrixCom

extern double Frobenius_s3(double *qqT, float *D);
extern double Frobenius_s3(double *qqT, double *D);
extern void p3update(double *D, double *W, double epsilon);
extern void p3update(double *D, double *invsqrtD, double *sqrtD, double *W, double epsilon);
extern void s3multi(double *A, double *B, double *AB);
extern void p3invsqrt(double *A, double *invsqrtA, double *sqrtA);
extern void p3invsqrt(double *A, double *invsqrtA);
extern int p3invsqrt(float *A, float *invsqrtA);
extern double p3RiemannianDistance(double *D, double *invsqrtF);
extern double p3RiemannianDistance(double *L);
extern float p3RiemannianDistance(float *D, float *invsqrtF);
extern void s3ABA(double *A, double *B, double *ABA);
extern void s3ULUT(double *L, double *UT, double *ULUT);
extern void s3eigenvec(double *A, double *evalue, double *evector);
extern void getcomplement1(double *a, double *b, double *c);
extern void getcomplement2 (double *U, double *V, double *W);
extern void s3inv(double *A, double *invA);
extern int s3inv(float *A, float *invA);
extern void s3multi(double *A, double *B, double *AB);
extern void differentiate_distance(double *D, double *F, double *L, double *dLdD);

extern void subtractRowAvg(float *X, int N, int P, float *X0);

extern void crossProduct(double *a, double *b, double *c);
extern void crossProduct(float *a, float *b, float *c);

extern void copyVector(float *v1, float *v2, int n);
extern void subtractVector(float *v1, float *v2, int n);
void normalizeVector(float *x, int n, double *norm);

extern int centerMatrixRow(float *X, int N, int P, float *avg);
extern int centerMatrixRow(float *X, int N, int P);

// Returns 1 if an error condition occurs, 0 otherwise
extern int avgRow(float *X, int N, int P, float *avg, char *rowmask);
extern int avgRow(double *X, int N, int P, double *avg, char *rowmask);
extern int avgRow(float *X, int N, int P, float *avg);
extern int avgRow(double *X, int N, int P, double *avg);
extern int varRow(float *X, int N, int P, float *avg, float *var);
extern int varRow(double *X, int N, int P, double *avg, double *var);
extern int varRow(double *X, int N, int P, double *avg, double *var, char *rowmask);
extern int varRow(float *X, int N, int P, float *avg, float *var, char *rowmask);
extern int ssdRow(float *X, int N, int P, float *avg, float *ssd, char *rowmask);
extern int ssdRow(float *X, int N, int P, float *avg, float *ssd);

// compute the Euclidian distance between two vectors r0 and r1
extern double euclideandistance(float *r0, float *r1, int n);
extern double xtAx(float *A, double *x, int p);
extern double vectorNorm(float *x, int n);
extern void normalizeVector(float *x, int n);
extern void transpose_matrix(float *A, int N,  int M);

extern float normalize(float *s, int n);
extern double normalize(double *s, int n);

extern int ComputeRank(double *M);
extern void s3eigenval(double *A, double *L);
extern double s3tr(double *A, double *B);
extern float s3tr(float *A, float *B);
extern void s3vec_to_mat(float *M, float *V);
extern void s3vec_to_mat(double *M, double *V);
extern void s3mat_to_vec(float *M, float *V);
extern void s3mat_to_vec(double *M, double *V);
extern void s3adjoint(double *A, double *ADJ);
extern void s3adjoint(float *A, float *ADJ);
extern double s3det(double *A);
extern float s3det(float *A);
void ds3det(double *A, double *B);
void ds3det(float *A, float *B);
extern double det3(double *A);
extern float det3(float *A);
extern float det4(float *A);
extern double det4(double *A);
extern double *inv3(double *A);
extern float *inv3(float *A);
extern void inv3(float *A, float *invA);
extern float *inv2(float *A);
extern double *inv2(double *A);
extern float *inv4(float *A);
extern double *inv4(double *A);
extern void multi(float *A,int iA,int jA,float *B,int iB,int jB,float *C);
extern void multi(double *A,int iA,int jA,double *B,int iB,int jB,double *C);
extern void multi(float *A,int iA,int jA, double *B,int iB,int jB,double *C);
extern void multi(double *A,int iA,int jA, float *B,int iB,int jB,float *C);
#endif

#ifndef _reslice

#define LIN 1
#define NEARN 2
#define SINC 3
#define CUBICSPLINE 4	

extern short *resliceImage(short *im1, DIM dim1, DIM dim2, float *T, int interpolation_method);

extern float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, float *xjit, float *yjit);

extern float partial_var(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np, float mu);

// you must initialize drand48 before using this function
extern unsigned char PNN(float x, float y, float z, unsigned char *array, int nx, int ny, int nz);

extern unsigned char nearestNeighbor(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np);
extern float nearestNeighbor(float x, float y, float z, float *array, int nx, int ny, int nz, int np);
extern short nearestNeighbor(float x, float y, float z, short *array, int nx, int ny, int nz, int np);


extern char *resliceImage(char *obj, int Onx, int Ony, float Odx, float Ody, int Tnx, int Tny, float Tdx, float Tdy, float *T);
extern short *resliceImage(short *im1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2, float *T);
extern short *resliceImage(float *im1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2, float *T);

extern short *resliceImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, int interpolation_method);

extern unsigned char *resliceImage(unsigned char *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T);

extern float *resliceImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *T, float *w);

extern short *resliceImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

extern unsigned char linearInterpolator(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np);
extern short linearInterpolator(float x, float y, float z, short *array, int nx, int ny, int nz, int np);
extern float linearInterpolator(float x, float y, float z, float *array, int nx, int ny, int nz, int np);
extern float linearInterpolator(float x, float y, float z, float *array, int nx, int ny, int nz, int np, float *w);
extern unsigned char linearInterpolator(float x, float y, float z, unsigned char *array, int nx, int ny, int nz, int np, float *w);


extern short *computeReslicedImage(short *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

extern float *computeReslicedImage(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp);

extern short *computeReslicedImage(float *im1, int nx1, int ny1, float dx1, float dy1,
int nx2, int ny2, float dx2, float dy2, float *Xwarp, float *Ywarp);

extern short *computeReslicedImage(short *im1, int nx1, int ny1, float dx1, float dy1,
int nx2, int ny2, float dx2, float dy2, float *Xwarp, float *Ywarp);


#endif

#ifndef _resize
extern short *resizeXYZ(short *image1,  DIM dim1, DIM dim2);

extern short *resizeXYZ(short *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

extern unsigned char *resizeXYZ(unsigned char *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

extern short *resizeXYZ(char *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

float *resizeXYZ(float *image1,
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

extern float *resizeXY(float *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);
extern short *resizeXY(short *image1, int nx1, int ny1, float dx1, float dy1, int nx2, int ny2, float dx2, float dy2);
#endif

#ifndef _getoption

extern int optind;
extern char *optarg;

struct option
{
	const char *name;
	int has_arg;
	int val;
};

extern int getoption(int argc, char **argv, struct option *options);
#endif

#ifndef _analyzeio
extern int extension_is_hdr(const char *filename);
extern void read_analyze_image(const char *filename, short *im);
extern char *read_analyze_image(const char *filename, DIM *dim, int *type, int v);
extern float read_dx(const char *hdrfile);
extern float read_dy(const char *hdrfile);
extern float read_dz(const char *hdrfile);
extern int read_nt(const char *hdrfile);
extern int read_nx(const char *hdrfile);
extern int read_ny(const char *hdrfile);
extern int read_nz(const char *hdrfile);
extern int read_datatype(char *hdrfile);
extern char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);
extern char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, int *type, int v);
extern char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type, int v);
extern char *read_analyze_image(const char *filename, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, int *type);
extern char *read_image(char *file,int n);
extern void get_analyze_file_names(const char *filename, char *basename_hdr, char *basename_img);
extern void read_analyze_hdr(struct dsr *hdr, char *filename);
extern void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz, short *dataType);
extern void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType);
extern void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz, short *dataType, int v);
extern void setDimensions(struct dsr hdr, int *nx, int *ny, int *nz, int *nt, float *dx, float *dy, float *dz, short *dataType, int v);
extern void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int dt, float dx, float dy, float dz);
extern void create_analyze_hdr(struct dsr *hdr, int nx, int ny, int nz, int nt, int datatype, float dx, float dy, float dz);
extern void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz); 
extern void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz); 
extern void write_analyze_image(const char *filename, unsigned char *im, int nx, int ny, int nz, float dx, float dy, float dz); 
extern void write_analyze_image(const char *filename, unsigned char *im, int nx, int ny, int nz, float dx, float dy, float dz, int v); 
extern void write_analyze_image(const char *filename, short *im, int nx, int ny, int nz, float dx, float dy, float dz,int v); 
extern void write_analyze_image(const char *filename, float *im, int nx, int ny, int nz, float dx, float dy, float dz,int v); 
#endif

#ifndef _swap
extern int bigEndian();
extern void swapByteOrder(char *in, int N);
extern void swapN(char *in, int N);
extern void swap_float_array( float *x, int n);
extern void swap_int_array( int *x, int n);
#endif

#ifndef _medianfilter
extern void medianFilter(float *image1, int nx, int ny, int nz, int Wx, int Wy, int Wz);
#endif

#ifndef _fileinfo
extern int isregular(const char *file);
extern int getFileSize(const char *file);
extern int checkWriteAccess(char *file);
extern int checkReadAccess(char *file);
extern int checkFileExistence(char *file);
extern int checkFileReadOK(char *file);
extern int checkFileWriteOK(char *file);
extern int check_F_R_permission(char *file);
#endif

#ifndef _histogram
extern int otsu(double *histogram, int numberOfBins);
extern double *findHistogram(short *im, int nv, int *nb, short &min, short &max);
extern double *findHistogram(short *im, int nv, int nb, int low, int high, int *bw_return);
extern double *findHistogram(short *im1, short *im2, int nv, int nb1, int nb2, int *bw1_r, int *bw2_r, int low1, int high1, int low2, int high2);
extern void trimExtremes(short *image, short *msk, int nv, float percent);
#endif

#ifndef _matrixops
extern void zeroVector(float *v, int n);
extern void zeroVector(char *v, int n);
extern void oneVector(float *v, int n);
extern void oneVector(char *v, int n);
extern void svd(float *At, int M, int N, float *Ut, float *V, float *S);
extern int zeroRowCol(float *A, int N, int n);
extern int setRowCol(float *A, int N, int n, float *a);
extern void projectVector(double *x, double *xpar, double *xper, float *Pz, int n);
extern float *projectionMatrix(double *X, int N, int p, int *rank);
extern float *projectionMatrix(float *X, int N, int p);
extern void mat_mat_trans(float *A,int Ar,int Ac,float *B,int Br, float *C);
extern float *diagATA_float(float *ATA, int n, char uplo);
extern float *AAT_float(float *A,int nr,int nc, char uplo);
extern void mat_trans_mat(float *A, int Ar, int Ac, float *B, int Bc, float *C);
extern void mat_trans_mat(double *A, int Ar, int Ac, double *B, int Bc, double *C);
#endif

#ifndef _nifti
extern int not_magical_nifti(const char *imagefilename);
extern char *read_nifti_image(const char *filename, nifti_1_header *hdr);
extern int same_nifti_image_size(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);
extern void read_nifti_image(const char *filename, unsigned char **im, nifti_1_header *hdr);
extern void read_nifti_image(const char *filename, short **im, nifti_1_header *hdr);
extern void print_NIFTI_hdr(const char *filename);
extern void print_NIFTI_hdr(nifti_1_header hdr);
extern nifti_1_header read_NIFTI_hdr(const char *filename);
extern int read_NIFTI_hdr(const char *filename, nifti_1_header *hdr);
extern nifti_1_header read_NIFTI_hdr(const char *filename, nifti1_extender *extender, char **extension);
extern void save_nifti_image(const char *filename, unsigned char *im, nifti_1_header *hdr);
extern void save_nifti_image(const char *filename, short *im, nifti_1_header *hdr);
extern void save_nifti_image(const char *filename, float *im, nifti_1_header *hdr);
extern void save_nifti_image(const char *filename, char *im, nifti_1_header *hdr);

// returns the orientations vectors xvec, yvec, and zvec in NIFTI's RAS system
extern void readOrientationVectorsFromFile(const char *filename, float *xvec, float *yvec, float *zvec);

extern void niftiFilename(char *filename, const char *path);
extern void swapniftiheader(nifti_1_header *hdr);
extern short *readNiftiImage(const char *filename, DIM *dim, int flg);
#endif


#ifndef _utils

// Set the nxn matrix A equal to the identity matrix
extern void set_to_I( float *A, int n);

extern void sobel_edge_x(short *in, float *out, int nx, int ny);
extern void sobel_edge_y(short *in, float *out, int nx, int ny);
extern void sobel_edge(short *in, float *out, int nx, int ny);
extern void sobel_edge(short *in, short *out, int nx, int ny);

extern int ccsize(short *im, int nv);

extern void copyarray(float *source, float *destination, int size);
extern void copyarray(short *source, char *destination, int size);
extern void copyarray(char *source, short *destination, int size);

extern void zeroarray(float *y, int size);
extern float diceindex(short *setA, short *setB, int n);
extern void remove_space(char *inp);
extern void orientationCodeConverter(int integerCode, char *stringCode);
extern int orientationCodeOK(char *stringCode);

extern void ijk2xyz(float *T, DIM dim);
extern void ijk2xyz(float *T, int nx, int ny, int nz, float dx, float dy, float dz);
extern void xyz2ijk(float *T, int nx, int ny, int nz, float dx, float dy, float dz);
extern void xyz2ijk(float *T, DIM dim);
extern void saveMatrix(float *A, int n, int m, char *filename);
extern void saveMatrix(short *A, int n, int m, char *filename);
extern float *readMatrix(int *n, int *m, char *filename);

extern float *readDataMatrix(char **imageList, int n, int p, short *mask);
extern float *readDataMatrix_nifti(char **imageList, int n, int p, short *mask);

extern short *readDataMatrixShort(char **imageList, int n, int p, short *mask);
extern short *readDataMatrixShort_nifti(char **imageList, int n, int p, short *mask);

extern short *readMask(const char *filename, int *nx, int *ny, int *nz);
extern short *readMask_nifti(const char *filename, int *nx, int *ny, int *nz);


extern void checkDimension(int N, char **imagefile, int nx, int ny, int nz);
extern void checkDimension_nifti(int N, char **imagefile, int nx, int ny, int nz);
extern void checkDimension(int N, char **imagefile, int *nx, int *ny, int *nz, float *dx, float *dy, float *dz);


extern void affineLSE(short *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
extern void affineLSE(char *msk, int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T);
extern float *affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp);
extern void affineLSE(char *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T);

extern void affineLSE(short *msk, int nx, int ny, float dx, float dy, float *Xwarp, float *Ywarp, float *T);

extern void extractArray(short *im, int nx, int ny, int nx0, int ny0, int Lx, int Ly, float *array);
extern void extractArray(short *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(float *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(unsigned char *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(short *im, int nx, int i, int j, int L, float *array);
extern void extractArray(short *im,int nx,int ny,int nz,int x0,int y0,int z0,short *x,short *y,short *z,int n,float *array);
extern int extractArray(short *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, short *array);
extern void extractArray(float *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(short *im, int nx, int ny, int nz, int np, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, float *array);
extern void extractArray(short *im, int nx, int ny, int nz, int nx0, int ny0, int nz0, int Lx, int Ly, int Lz, short *array);

// This function extracts the "filename" from the full "path" string.
// For example, if path="/home/babak/testimages/test.img", then filname="test.img".
extern void getfilename(char *filename, const char *path);

extern void printMatrix(int *mat, int n, int p, const char *s, FILE *fp);
extern void printMatrix(float *mat, int n, int p, const char *s, FILE *fp);
extern void printMatrix(double *mat, int n, int p, const char *s, FILE *fp);
extern void get_temp_filename(char *filename);
extern void mask_and_save(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
extern void mask_and_save_nii(const char *inputfile, const char *outputfile, short *mask, short *masked_image, int nbv, float FWHM);
extern void read_transpose_save(char *inputfile, char *outputfile, int nr, int v);
extern void centerOfMass(short *im, int nx, int ny, int nz, float dx, float dy, float dz, float *CM);
#endif

extern float *smoothY(float *image, int nx, int ny, int nz, float sd);
extern float *smoothZ(float *image, int nx, int ny, int nz, float sd);

#endif
