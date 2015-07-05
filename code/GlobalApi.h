
#include "DifImage.h"
#include "GlobalParam.h"

#define NORM_GRAY 255

#define MIN_DOUBLE 1.0e-15

#define TRANSFORM

#define CBASE 1.5
#define LOGC 0.40547
#define PI 3.1415926535
#define C1 0
#define C2 0

#define W_SGN 1.2

#define ALPHA_BV 0.01 // coordinate for BV
#define BETA_INV 1.0 // coordinate for invBase
#define DISCTYPT_TV 2 // 0 naive; 1 sample; 2 sophi 
#define TP 4.0 // double
#define TASK 1 // 0: denoise; 1: inpaint

// LPdeO.cpp
void train (string inFilePath, string outFilePath, string geneFilePath, string inMaskPath, string taskType, double *a, GlobalParam_T &globalParam);
void test (string testFilePath, string geneFilePath, string testMaskPath, double *a, GlobalParam_T &globalParam);
void showStepResult(double *val);
void initA(double *a, double aLength);

// Optimization.cpp
void minimize (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam);
void grad (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *g, GlobalParam_T &globalParam);
double linMin (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *d, double *da, GlobalParam_T &globalParam);
double objVal (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *d, double w, GlobalParam_T &globalParam);

// PdeSolver.cpp
void pdeSolver (DifImage_T &inImage, DifImage_T &maskImage, double *pGenImageData1, double *a, double dt);
void pdeSolver (DifImage_T &inImage, DifImage_T &maskImage, DifImage_T *pGen, double *a, double dt);
void pdeStep (DifImage_T &gen, DifImage_T &nextGen, DifImage_T &inImage, DifImage_T &maskImage, double *a, double dt);
void setBoundary (double *image, double value, int width, int height);

// AdjointPdeSolver.cpp
void pdeAdjSolver (DifImage_T &inImage, DifImage_T &maskImage, DifImage_T &outImage, DifImage_T *pGen, double *pLagMul, double *a, double dt);
void pdeAdjStep (DifImage_T &inImage, DifImage_T &maskImage, DifImage_T &gen, double *pLagMul, double *pPrevLagMul, double *a, double dt);
void geneSigma (DifImage_T &O, DifImage_T &maskImage, double *sigma, double *a, int p, int q);
void geneSigmaBase (DifImage_T &O, double *sigmaBase, bool *sigmaBaseNonzero, int p, int q);
void geneSigmaTV(DifImage_T &O, double *sigmaTV, bool &sigmaTVNonzero, int p, int q);

// MathOperation.cpp
void multiply (double *A, double *b, double *Ab, int n);
double innerProduct (double *a, double *b, int n);

// Invariants.cpp
void geneInv(DifImage_T &O, DifImage_T &inImage, double *invBase);
void geneBV(DifImage_T &O, DifImage_T &inImage, DifImage_T &maskImage, double *bv);
void geneSgn(DifImage_T &O, DifImage_T &inImage, double *sgn);
void geneB(double *image, double *pB, int idx, int width, int height);
void geneLap(DifImage_T &O, double *lap);
void geneLambda(DifImage_T &O, double *lambda);

// DifImage.cpp
void difImage (double *image, double *difImage, int width, int height);
void difImage (double *image, double *difImage, int i, int j, int width, int height);

// FileIO.cpp
bool imageRead(DifImage_T &inImage, string imageFile, int gapWidth);
bool imageWrite(double *val, string imageFile, int width, int height, int gapWidth);

// Initial.cpp
void init (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam);
void stepInit (DifImage_T *pGen, DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, int t, GlobalParam_T &globalParam);
void qcqp (double *F, double *b, double c, double *a);
double lineFunW (double *F, double *b, double c, double *a, double *d);
double funW (double *F, double *b, double c, double *a, double *d, double w);
void restart (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam);

///////////////////////////////////////////////////////////////////
inline void setZeroValue (double *v, int n)
{
	for (int i = 0; i < n; i++)
	{
		v[i] = 0;
	}
}

inline void setZeroBoundary (double *image, int width, int height)
{
	double *pI, *pIEnd;
	pI = image;
	pIEnd = pI + width;
	while (pI < pIEnd)
	{
		*pI = 0;
		pI++;
	}
	pI = image + (height - 1) * width;
	pIEnd = pI + width;
	while (pI < pIEnd)
	{
		*pI = 0;
		pI++;
	}
	pI = image;
	pIEnd = pI + (height - 1) * width;
	while (pI < pIEnd)
	{
		pI += (width - 1);
		*pI = 0;
		pI++;
		*pI = 0;
	}
	return;
}

//storage of coefficient matrices
//i   j  r
//0   0  0
//0   1  1
//1   0  2
//0   2  3
//1   1  4
//2   0  5
inline int ind (int i, int j)      //index mapping of f_(i, j)
{
	int index;
	if (i == 0 && j == 0)
	{
		index = 0;
	}
	else if (i == 0 && j == 1)
	{
		index = 1;
	}
	else if (i == 1 && j == 0)
	{
		index = 2;
	}
	else if (i == 0 && j == 2)
	{
		index = 3;
	}
	else if (i == 1 && j == 1)
	{
		index = 4;
	}
	else if (i == 2 && j == 0)
	{
		index = 5;
	}
	else
	{
		index = -1;
	}

	return index;
}

inline void invInd (int index, int &i, int &j)      //inverse index mapping of f_(i, j)
{
	switch (index)	{
	case 0: i = 0; j = 0; break;
	case 1: i = 0; j = 1; break;
	case 2:	i = 1; j = 0; break;
	case 3: i = 0; j = 2; break;
	case 4: i = 1; j = 1; break;
	case 5: i = 2; j = 0; break;
	default: i = -1; j = -1;
	}
	return;
}

/*inline int sgn (double x)
{
	int sign = 0;
	if (x > 0)
	{
		sign = 1;
	}
	else if (x < 0)
	{
		sign = -1;
	}
	return sign;
}*/