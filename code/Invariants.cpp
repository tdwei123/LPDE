
#include "StdAfx.h"
#include "GlobalApi.h"
#include "DifImage.h"

//************************************************************************
//*  Function Name: geneInv
//*  Function Description: 
//*      computes invariant bases inv_i(t).
//*  Arguments: 
//*      [IN] : DifImage_T &O
//*      [OUT] : double *invBase
//* 
//*  Return Value: void 
//*      none.
//* 
//*  Last Modified on: 2009-8-3 13:46:52 by Risheng Liu
//************************************************************************

void geneInv(DifImage_T &O, DifImage_T &inImage, double *invBase)
{
	int j, jStart;
	int n = O.getWidth() * O.getHeight();

	double *O0 = O.getVal(0, 0);
	double *Ox = O.getVal(1, 0);
	double *Oy = O.getVal(0, 1);
	double *Oxx = O.getVal(2, 0);
	double *Oxy = O.getVal(1, 1);
	double *Oyy = O.getVal(0, 2);

	double *f0 = inImage.getVal(0, 0); // modified by Risheng 2009-11-12
	double *lap = new double[n];

	geneLap(O, lap);

	// inv_0
	for (j = 0; j < n; j++)
	{
		if (TASK == 0)
		{
			invBase[j] = f0[j];
			//invBase[j] = 1;
		}
		else if (TASK == 1)
		{
			invBase[j] = 0;
		}
	}

	// inv_1
	jStart = n;
	for (j = 0; j < n; j++)
	{
		invBase[jStart + j] = O0[j];
	}

	// inv_2
	jStart = jStart + n;
	for (j = 0; j < n; j++)
	{
		invBase[jStart + j] = Ox[j] * Ox[j] + Oy[j] * Oy[j];
	}

	// inv_3
	jStart = jStart + n;
	for (j = 0; j < n; j++)
	{
		//invBase[jStart + j] = Oxx[j] + Oyy[j];
		invBase[jStart + j] = lap[j];
	}

	// inv_4
	jStart = jStart + n;
	for (j = 0; j < n; j++)
	{
		invBase[jStart + j] = (Ox[j] * Oxx[j] + Oy[j] * Oxy[j]) * Ox[j] 
			+ (Ox[j] * Oxy[j] + Oy[j] * Oyy[j]) * Oy[j];
	}

	// inv_5
	jStart = jStart + n;
	for (j = 0; j < n; j++)
	{
		invBase[jStart + j] = Oxx[j] * Oxx[j] + 2 * Oxy[j] * Oxy[j] + Oyy[j] * Oyy[j];
	}

	delete []lap;

}

void geneBV(DifImage_T &O, DifImage_T &inImage, DifImage_T &maskImage, double *bv)
{
	int j, jStart;
    int width = O.getWidth();
	int height = O.getHeight();
	int n = width * height;
	double eps = 1e-6;

	// for case 0
	double tmp1, tmp2;

	double *O0 = O.getVal(0, 0);			
		
	double *Ox = O.getVal(1, 0);
	double *Oy = O.getVal(0, 1);
	double *Oxx = O.getVal(2, 0);
	double *Oxy = O.getVal(1, 1);
	double *Oyy = O.getVal(0, 2);

	double *mask = maskImage.getVal(0, 0);

	// for case 1
	int x, y, yBase, yForward, yBackward;
	
	double *b1 = new double[n];//b1 = bi-1/2j
	double *b2 = new double[n];//b2 = bi+1/2j
	double *b3 = new double[n];//b3 = bij+1/2
	double *b4 = new double[n];//b4 = bij-1/2 
	
	geneB(O0, b1, 1, width, height);
	geneB(O0, b2, 2, width, height);
	geneB(O0, b3, 3, width, height);
	geneB(O0, b4, 4, width, height);

	// for case 2
	double *b5 = new double[n];//b5 = bi-1/2j-1/2
	double *b6 = new double[n];//b6 = bi-1/2j+1/2
	double *b7 = new double[n];//b7 = bi+1/2j-1/2
	double *b8 = new double[n];//b8 = b1+1/2j+1/2

	geneB(O0, b5, 5, width, height);
	geneB(O0, b6, 6, width, height);
	geneB(O0, b7, 7, width, height);
	geneB(O0, b8, 8, width, height);

	double *lambda = new double[n];//lambda
	geneLambda(O, lambda);

	double *sgn = new double[n];
	geneSgn(O, inImage, sgn);

	switch (DISCTYPT_TV)
	{
	case 0:
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] * Oyy[j] + Oy[j] * Oy[j] * Oxx[j] - 2 * Ox[j] * Oy[j] * Oxy[j];
			tmp2 = pow((Ox[j] * Ox[j] + Oy[j] * Oy[j]), 1.5) + eps;
			bv[j] = tmp1 / tmp2;
		}
		setZeroBoundary(bv, width, height);
		break;

	case 1:
		yBase = 0;
		yForward = width;
		for (y = 1; y < height; y++)
		{
			yBackward = yBase;
			yBase = yForward;
			yForward += width;
			for (x = 1; x < width - 1; x++)
			{
				bv[x + yBase] = b1[x + yBase] * O0[x - 1 + yBase] + b2[x + yBase] * O0[x + 1 + yBase]
					+ b3[x + yBase] * O0[x + yForward] + b4[x + yBase] * O0[x + yBackward]
					- (b1[x + yBase] + b2[x + yBase] + b3[x + yBase] + b4[x + yBase]) * O0[x + yBase];
			}
		}
		setZeroBoundary(bv, width, height);
		break;

	case 2:
		yBase = 0;
		yForward = width;
		for (y = 1; y < height; y++)
		{
			yBackward = yBase;
			yBase = yForward;
			yForward += width;
			for (x = 1; x < width - 1; x++)
			{
				tmp1 = b1[x + yBase] * O0[x - 1 + yBase] + b2[x + yBase] * O0[x + 1 + yBase]
					+ b3[x + yBase] * O0[x + yForward] + b4[x + yBase] * O0[x + yBackward]
					- (b1[x + yBase] + b2[x + yBase] + b3[x + yBase] + b4[x + yBase]) * O0[x + yBase];
				tmp2 = (b5[x + yBase] * O0[x - 1 + yBackward] + b6[x + yBase] * O0[x - 1 + yForward]
					+ b7[x + yBase] * O0[x + 1 + yBackward] + b8[x + yBase] * O0[x + 1 + yForward]
					- (b5[x + yBase] + b6[x + yBase] + b7[x + yBase] + b8[x + yBase]) * O0[x + yBase]) / 2;
				//bv[x + yBase] = (tmp1 + tmp2) / 2;
				bv[x + yBase] = lambda[x + yBase] * tmp1 + (1 - lambda[x + yBase]) * tmp2;
			}
		}
		setZeroBoundary(bv, width, height);
		break;
	default:
		cout<<"Error TV Type!"<<endl;
	}

	jStart = n;
	for (j = 0; j < n; j++)
	{
		bv[jStart + j] = W_SGN * sgn[j];
		/*if (TASK == 0)
		{
			bv[jStart + j] = W_SGN * sgn[j];//weight
		}
		else if (TASK == 1)
		{
			if (mask[j] == 0)
			{
				bv[jStart + j] = 1;
			}
			else 	
			{
				bv[jStart + j] = W_SGN * sgn[j];//weight
			}
		}*/
	}
	
	delete []b1;
	delete []b2;
	delete []b3;
	delete []b4;
	delete []b5;
	delete []b6;
	delete []b7;
	delete []b8;

	delete []lambda;
	delete []sgn;

	return;
}

void geneLambda(DifImage_T &O, double *lambda)
{
	int j;
	int width = O.getWidth();
	int height = O.getHeight();
	int n = width * height;

	double *O0 = O.getVal(0, 0);
	double *Ox = O.getVal(1, 0);
	double *Oy = O.getVal(0, 1);
	double theta;

	double eps = 1e-7;

	for (j = 1; j < n; j++)
	{
		theta = fabs(atan((Ox[j])/(Oy[j] + eps)));
		if (theta > 0 && theta < PI / 4)
		{
			lambda[j] = 1 - 4 * theta / PI;
		}
		else
		{
			lambda[j] = 4 * (theta - PI / 4) / PI; 
		}
		
	}
	setZeroBoundary(lambda, width, height);
}

//b1 = bi-1/2j; b2 = bi+1/2j; b3 = bij+1/2; b4 = bij-1/2;
//b5 = bi-1/2j-1/2; b6 = bi-1/2j+1/2; b7 = bi+1/2j-1/2; b8 = b1+1/2j+1/2;
void geneB(double *image, double *pB, int idx, int width, int height)
{
	int x, y, yBase, yForward, yBackward;
	int n = width * height;
	double tmp1;
	double tmp2;
	double eps = 1e-6;

	yBase = 0;
	yForward = width;

	if (idx == 1)//b1 = bi-1/2j
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + yBase] - image[x - 1 + yBase]) * (image[x + yBase] - image[x - 1 + yBase]);
				tmp2 = (image[x + yForward] + image[x - 1 + yForward] - image[x + yBackward] - image[x - 1 + yBackward])
					* (image[x + yForward] + image[x - 1 + yForward] - image[x + yBackward] - image[x - 1 + yBackward]);
				pB[x + yBase] = 1 / (sqrt(tmp1 + (tmp2 / 16)) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 2)//b2 = bi+1/2j
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + 1 + yBase] - image[x + yBase]) * (image[x + 1 + yBase] - image[x + yBase]);
				tmp2 = (image[x + 1 + yForward] + image[x + yForward] - image[x + 1 + yBackward] - image[x + yBackward])
					* (image[x + 1 + yForward] + image[x + yForward] - image[x + 1 + yBackward] - image[x + yBackward]);
				pB[x + yBase] = 1 / (sqrt(tmp1 + (tmp2 / 16)) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 3)//b3 = bij+1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + yForward] - image[x + yBase]) * (image[x + yForward] - image[x + yBase]);
				tmp2 = (image[x + 1 + yForward] + image[x + 1 + yBase] - image[x - 1 + yForward] - image[x - 1 + yBase])
					* (image[x + 1 + yForward] + image[x + 1 + yBase] - image[x - 1 + yForward] - image[x - 1 + yBase]);
				pB[x + yBase] = 1 / (sqrt(tmp1 + (tmp2 / 16)) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 4)//b4 = bij-1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + yBase] - image[x + yBackward]) * (image[x + yBase] - image[x + yBackward]);
				tmp2 = (image[x + 1 + yBase] + image[x + 1 + yBackward] - image[x - 1 + yBase] - image[x - 1 + yBackward])
					* (image[x + 1 + yBase] + image[x + 1 + yBackward] - image[x - 1 + yBase] - image[x - 1 + yBackward]);
				pB[x + yBase] = 1 / (sqrt(tmp1 + (tmp2 / 16)) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 5)//b5 = bi-1/2j-1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x - 1 + yBase] - image[x + yBackward]) * (image[x - 1 + yBase] - image[x + yBackward]);
				tmp2 = (image[x + yBase] - image[x - 1 + yBackward]) * (image[x + yBase] - image[x - 1 + yBackward]);
				pB[x + yBase] = 1 / (sqrt((tmp1 / 2) + (tmp2 / 2) ) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 6)//b6 = bi-1/2j+1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x - 1 + yForward] - image[x + yBase]) * (image[x - 1 + yForward] - image[x + yBase]);
				tmp2 = (image[x + yForward] - image[x - 1 + yBase]) * (image[x + yForward] - image[x - 1 + yBase]);
				pB[x + yBase] = 1 / (sqrt((tmp1 / 2) + (tmp2 / 2) ) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 7)//b7 = bi+1/2j-1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + yBase] - image[x + 1 + yBackward]) * (image[x + yBase] - image[x + 1 + yBackward]);
				tmp2 = (image[x + 1 + yBase] - image[x + yBackward]) * (image[x + 1 + yBase] - image[x + yBackward]);
				pB[x + yBase] = 1 / (sqrt((tmp1 / 2) + (tmp2 / 2) ) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	else if (idx == 8)//b8 = b1+1/2j+1/2
	{
		for (y = 1; y < height - 1; y++)
		{
			yBackward = yBase;
		    yBase = yForward;
		    yForward += width;
			for(x  = 1; x < width - 1; x++)
			{
				tmp1 = (image[x + 1 + yForward] - image[x + yBase]) * (image[x + 1 + yForward] - image[x + yBase]);
				tmp2 = (image[x + yForward] - image[x + 1 + yBase]) * (image[x + yForward] - image[x + 1 + yBase]);
				pB[x + yBase] = 1 / (sqrt((tmp1 / 2) + (tmp2 / 2)) + eps); 
			}
		}
		setZeroBoundary(pB, width, height);
	}
	
}

void geneSgn(DifImage_T &O, DifImage_T &inImage, double *sgn)
{
	int j;
	int width = O.getWidth();
	int height = O.getHeight();
	int n = width * height;
	double *O0 = O.getVal(0, 0);
	double *f0 = inImage.getVal(0, 0);
	double eps = 1e-6;

	for (j = 0; j < n; j++)
	{
		sgn[j] = (f0[j] - O0[j]) / (fabs(f0[j] - O0[j]) + eps);
	}
	setZeroBoundary(sgn, width, height);
}

void geneLap(DifImage_T &O, double *lap)
{
    int width = O.getWidth();
	int height = O.getHeight();
	int n = width * height;

	double *image = O.getVal(0, 0);
	//double *Oxx = O.getVal(2, 0);
	//double *Oyy = O.getVal(0, 2);

	int x, y, yBase, yForward, yBackward;

	yBase = 0;
	yForward = width;
		
	for (y = 1; y < height; y++)
	{
		yBackward = yBase;
		yBase = yForward;
		yForward += width;
		for (x = 1; x < width - 1; x++)
		{
			lap[x + yBase] = (image[x + yForward] + image[x + yBackward] + image[x + 1 + yBase] + image[x - 1 + yBase]
				+ image[x + 1 + yForward] + image[x + 1 + yBackward] + image[x - 1 + yForward] 
				+ image[x - 1 + yBackward] - 8 * image[x + yBase]) / 3;
			//lap[x + yBase] = Oxx[x + yBase] + Oyy[x + yBase];
		}	
	}
	setZeroBoundary(lap, width, height);
}