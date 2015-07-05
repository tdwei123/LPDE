
#include "StdAfx.h"
#include "GlobalApi.h"
#include "DifImage.h"

DifImage_T::DifImage_T (int width, int height)
{
	i_width = width;
	i_height = height;
	int n = i_width * i_height;
	if (n == 0)
	{
		i_pVal = NULL;
	}
	else
	{
		i_pVal = new double[6 * n];
	}
}

DifImage_T::~DifImage_T(void)
{
	if (i_pVal != NULL)
	{
		delete []i_pVal;
	}
	i_width = 0;
	i_height = 0;
}

void DifImage_T::setArea (int width, int height)
{
	if (i_width != width || i_height != height)
	{
		i_width = width;
		i_height = height;
		if (i_pVal != NULL)
		{
			delete []i_pVal;
		}
		i_pVal = new double[6 * i_width * i_height];
	}
}

void DifImage_T::setVal (double *pImageVal)
{
	int n = i_width * i_height;
	if (n != 0)
	{
		difImage(pImageVal, i_pVal, i_width, i_height);
	}
	return;
}

void DifImage_T::setVal (DifImage_T &img)
{
	int n = i_width * i_height;
	if (i_width != img.getWidth() || i_height != img.getHeight())
	{
		i_width = img.getWidth();
		i_height = img.getHeight();
		if (i_pVal != NULL)
		{
			delete []i_pVal;
		}
		n = i_width * i_height;
		i_pVal = new double[6 * n];
	}
	double *pVal = img.getVal(0, 0);
	for (int i = 0; i < 6 * n; i++)
	{
		i_pVal[i] = pVal[i];
	}

	return;
}

double DifImage_T::getVal (int i, int j, int x, int y)
{
	int n = i_width * i_height;
	return i_pVal[ind(i, j) * n + x + y * i_width];
}

double DifImage_T::getVal (int i, int j, int pos)
{
	int n = i_width * i_height;
	return i_pVal[ind(i, j) * n + pos];
}

double *DifImage_T::getVal (int i, int j)
{
	int n = i_width * i_height;
	return i_pVal + ind(i, j) * n;
}

//************************************************************************
//*  Function Name: difImage 
//*  Function Description: 
//*      calculate the image differential till second order 
//*  Arguments: 
//*      [IN] : double *image
//*      [OUT] : double *difImage (the length of difImage is 6 * n)
//*      [IN] : int width
//*      [IN] : int height
//* 
//*  Return Value: void 
//*       none.
//* 
//*  Last Modified on: 2009-8-3 10:23:36 by Risheng Liu
//************************************************************************

void difImage (double *image, double *difImage, int width, int height)
{
	int x, y, yBase, yForward, yBackward;
	int n = width * height;
	double *pDifImage;

	pDifImage = difImage;
	// i == 0 && j == 0
	for (x = 0; x < n; x++)
	{
		pDifImage[x] = image[x];
	}
	setZeroBoundary(pDifImage, width, height);
	pDifImage += n;
	// dy: i == 0 && j == 1
	yBase = 0;
	yForward = width;
	for (y = 1; y < height - 1; y++)
	{
		yBackward = yBase;
		yBase = yForward;
		yForward += width;
		for (x = 1; x < width - 1; x++)
		{
			pDifImage[x + yBase] = (image[x + yForward] - image[x + yBackward]) / 2;
		}
	}
	setZeroBoundary(pDifImage, width, height);
	pDifImage += n;
	// dx: i == 1 && j == 0
	yBase = 0;
	for (y = 1; y < height - 1; y++)
	{
		yBase += width;
		for (x = 1; x < width - 1; x++)
		{
			pDifImage[x + yBase] = (image[x + 1 + yBase] - image[x - 1 + yBase]) / 2;
		}
	}
	setZeroBoundary(pDifImage, width, height);
	pDifImage += n;
	// dyy: i == 0 && j == 2
	yBase = 0;
	yForward = width;
	for (y = 1; y < height - 1; y++)
	{
		yBackward = yBase;
		yBase = yForward;
		yForward += width;
		for (x = 1; x < width - 1; x++)
		{
			pDifImage[x + yBase] = image[x + yBackward] + image[x + yForward] - image[x + yBase] - image[x + yBase];
		}
	}
	setZeroBoundary(pDifImage, width, height);
	pDifImage += n;
	// dxy i == 1 && j == 1
	yBase = 0;
	yForward = width;
	for (y = 1; y < height - 1; y++)
	{
		yBackward = yBase;
		yBase = yForward;
		yForward += width;
		for (x = 1; x < width - 1; x++)
		{
			pDifImage[x + yBase] = (image[x + 1 + yForward] + image[x - 1 + yBackward]
				- image[x - 1 + yForward] - image[x + 1 + yBackward]) / 4;
		}
	}
	setZeroBoundary(pDifImage, width, height);
	pDifImage += n;
	// dxx i == 2 && j == 0
	yBase = 0;
	for (y = 1; y < height - 1; y++)
	{
		yBase += width;
		for (x = 1; x < width - 1; x++)
		{
			pDifImage[x + yBase] = image[x + 1 + yBase]
				+ image[x - 1 + yBase] - image[x + yBase] - image[x + yBase];
		}
	}
	setZeroBoundary(pDifImage, width, height);
}

//************************************************************************
//*  Function Name: difImage 
//*  Function Description: 
//*      calculate the image differential, indexed by (i,j), for adjoint equation
//*  Arguments: 
//*      [IN] : double *image
//*      [IN] : double *difImage (the length of difImage is n)
//*      [IN] : int i
//*      [IN] : int j
//*      [IN] : int width
//*      [IN] : int height
//* 
//*  Return Value: void 
//*       none.
//* 
//*  Last Modified on: 2009-8-3 10:17:57 by Risheng Liu
//************************************************************************

void difImage (double *image, double *difImage, int i, int j, int width, int height)
{
	int x, y, yForward, yBase, yBackward;
	int n = width * height;
	if (i == 0 && j == 0)
	{
		for (x = 0; x < n; x++)
		{
			difImage[x] = image[x];
		}
	}
	else if (i == 0 && j == 1)
	{
		yBase = 0;
		yForward = width;
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
			yBase = yForward;
			yForward += width;
			for (x = 1; x < width - 1; x++)
			{
				difImage[x + yBase] = (image[x + yForward] - image[x + yBackward]) / 2;
			}
		}
	}
	else if (i == 1 && j == 0)
	{
		yBase = 0;
		for (y = 1; y < height - 1; y++)
		{
			yBase += width;
			for (x = 1; x < width - 1; x++)
			{
				difImage[x + yBase] = (image[x + 1 + yBase] - image[x - 1 + yBase]) / 2;
			}
		}
	}
	else if (i == 0 && j == 2)
	{
		yBase = 0;
		yForward = width;
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
			yBase = yForward;
			yForward += width;
			for (x = 1; x < width - 1; x++)
			{
				difImage[x + yBase] = image[x + yBackward] + image[x + yForward] - image[x + yBase] - image[x + yBase];
			}
		}
	}
	else if (i == 1 && j == 1)
	{
		yBase = 0;
		yForward = width;
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
			yBase = yForward;
			yForward += width;
			for (x = 1; x < width - 1; x++)
			{
				difImage[x + yBase] = (image[x + 1 + yForward] + image[x - 1 + yBackward]
					- image[x - 1 + yForward] - image[x + 1 + yBackward]) / 4;
			}
		}
	}
	else if (i == 2 && j == 0)
	{
		yBase = 0;
		for (y = 1; y < height - 1; y++)
		{
			yBase += width;
			for (x = 1; x < width - 1; x++)
			{
				difImage[x + yBase] = image[x + 1 + yBase]
					+ image[x - 1 + yBase] - image[x + yBase] - image[x + yBase];
			}
		}
	}
	setZeroBoundary(difImage, width, height);
}