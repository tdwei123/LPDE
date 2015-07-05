
#include "StdAfx.h"
#include "GlobalApi.h"
#include "DifImage.h"

void pdeSolver (DifImage_T &inImage, DifImage_T &maskImage, double *pGenImageData1, double *a, double dt)
{
	int i;
	int width = inImage.getWidth();
	int height = inImage.getHeight();
	int n = width * height;
	int totalTime = int(TP / dt);
	DifImage_T *pGen = new DifImage_T[totalTime + 1];
	for (i = 0; i <= totalTime; i++)
	{
		pGen[i].setArea(width, height);
	}
	pdeSolver(inImage, maskImage, pGen, a, dt);
	double *pGenImageData = pGen[totalTime].getVal(0, 0);
	for (i = 0; i < n; i++)
	{
		pGenImageData1[i] = pGenImageData[i];
	}
	delete []pGen;
}

void pdeSolver(DifImage_T &inImage, DifImage_T &maskImage, DifImage_T *pGen, double *a, double dt)
{
	int totalTime = int(TP / dt);
	int i;
	double *pInitImage = inImage.getVal(0, 0);
	pGen[0].setVal(pInitImage);
	for (i = 0; i < totalTime; i++)
	{
		pdeStep(pGen[i], pGen[i + 1], inImage, maskImage, &a[i * INV_NUM], dt);
	}
	return;
}

void pdeStep(DifImage_T &gen, DifImage_T &nextGen, DifImage_T &inImage, DifImage_T &maskImage, double *a, double dt)
{
	int i, k, idx;
	int width = gen.getWidth();
	int height = gen.getHeight();
	int n = width * height;
	double *inv = new double[INV_NUM * n];
	double *bv = new double[2 * n];
	double *L = new double[n];
	double *pCurImage;

	double *mask = maskImage.getVal(0, 0);

	geneInv(gen, inImage, inv);

	// dO/dt = Lo;
	pCurImage = gen.getVal(0, 0);
	for (k = 0; k < n; k++)
	{
		L[k] = 0;
	}

	idx = 0;
	for (i = 0; i < INV_NUM; i++)
	{
		for (k = 0; k < n; k++)
		{
			//////////////////////////////////////////
			if (TASK == 0)
			{
				L[k] = L[k] + a[i] * inv[idx];
				idx++;
			}
			else if (TASK == 1)
			{
				if (mask[k] == 0)
				{
					L[k] = L[k] + a[i] * inv[idx];
					idx++;
				} 
				else
				{
					idx++;
				}
			}
			//////////////////////////////////
		}
	}

    // generate BV
	geneBV(gen, inImage, maskImage, bv);
	for (k = 0; k < n; k++)
	{
		L[k] = BETA_INV * L[k] + ALPHA_BV * (bv[k] + bv[n + k]);// bv[k]--tv, bv[k + n]--sgn
	}

	for (k = 0; k < n; k++)
	{
		L[k] = pCurImage[k] + dt * L[k];
	}
	setBoundary(L, 0, width, height);
	nextGen.setVal(L);

	delete []inv;
	delete []L;
	delete []bv;
	return;
}

//************************************************************************
//*  Function Name: setBoundary 
//*  Function Description: 
//*      assign the boundary value, for LPDE, we need all the image with zero boundary as the boundary condition. 
//*  Arguments: 
//*      [IN] : double *image
//*      [IN] : double value
//*      [IN] : int width
//*      [IN] : int height
//* 
//*  Return Value: void 
//*      none.
//* 
//*  Last Modified on: 2009-8-5 0:08:21 by Risheng Liu
//************************************************************************

void setBoundary (double *image, double value, int width, int height)
{
	int k;
	for (k = 0; k < width; k++)
	{
		image[k] = value;
		image[k + (height - 1) * width] = value;
	}
	for (k = 1; k < height - 1; k++)
	{
		image[k * width] = value;
		image[width - 1 + k * width] = value;
	}
	return;
}