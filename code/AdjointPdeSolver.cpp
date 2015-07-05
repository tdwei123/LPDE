
#include "StdAfx.h"
#include "GlobalApi.h"
#include "DifImage.h"

void pdeAdjSolver(DifImage_T &inImage, DifImage_T &maskImage, DifImage_T &outImage, DifImage_T *pGen, double *pLagMul, double *a, double dt)
{
	int totalTime = int(TP / dt);
	int i, timeInd, iStart;
	int width = inImage.getWidth();
	int height = inImage.getHeight();
	int n = width * height;

	pdeSolver(inImage, maskImage, pGen, a, dt);
	double *pGen1 = pGen[totalTime].getVal(0, 0);
	double *pOut = outImage.getVal(0 ,0);

	timeInd = (totalTime - 1) * n;
	for (i = 0; i < n; i++) // phi(1) = O~ - O(1)
	{
		pLagMul[timeInd + i] = pOut[i] - pGen1[i];
	}
	setBoundary(&pLagMul[timeInd], 0, width, height);

	iStart = (totalTime - 1) * INV_NUM;
	// pdeAdjStep(inImage, pGen[totalTime], &pLagMul[timeInd], &pLagMul[timeInd - n], &a[iStart], dt);
	for (i = totalTime - 1; i > 0; i--)
	{
		pdeAdjStep(inImage, maskImage, pGen[i], &pLagMul[timeInd], &pLagMul[timeInd - n], &a[iStart], dt);
		iStart = iStart - INV_NUM;
		timeInd = timeInd - n;
	}
	return;
}

void pdeAdjStep(DifImage_T &inImage, DifImage_T &maskImage, DifImage_T &gen, 
				double *pLagMul, double *pPrevLagMul, double *a, double dt)
{
	int i, p, q, k;
	int width = inImage.getWidth();
	int height = inImage.getHeight();
	int n = width * height;
	double *sigma = new double[n];
	double *tmp = new double[n];
	double *d_tmp = new double[n];
	double *L = new double[n];

	setZeroValue(L, n);

	for (i = 0; i < 6; i++)
	{
		invInd(i, p, q);
		geneSigma(gen, maskImage, sigma, a, p, q);
		
		for (k = 0; k < n; k++)
		{
			tmp[k] = sigma[k] * pLagMul[k];
		}

		difImage(tmp, d_tmp, p, q, width, height);
		if ((p == 1 && q == 0) || (p == 0 && q == 1))
		{
			for (k = 0; k < n; k++)
			{
				L[k] = L[k] - d_tmp[k];
			}
		}
		else
		{
			for (k = 0; k < n; k++)
			{
				L[k] = L[k] + d_tmp[k];
			}
		}
	}
	for (k = 0; k < n; k++)
	{
		pPrevLagMul[k] = pLagMul[k] + dt * L[k];
	}
	setBoundary(pPrevLagMul, 0, width, height);

	delete []sigma;
	delete []tmp;
	delete []d_tmp;
	delete []L;

	return;
}	

void geneSigma(DifImage_T &O, DifImage_T &maskImage, double *sigma, double *a, int p, int q)
{
	int i, k, idx;
	int n = O.getWidth() * O.getHeight();
	double *inv = new double[INV_NUM * n];
	double *sigmaBase = new double[INV_NUM * n];
	bool sigmaBaseNonZero[INV_NUM];
	double *sigmaTV = new double[n];
	bool sigmaTVNonzero = true;

	double *mask = maskImage.getVal(0, 0);

	for (k = 0; k < n; k++)
	{
		sigma[k] = 0;
	}
	geneSigmaBase(O, sigmaBase, sigmaBaseNonZero, p, q);

	idx = 0;
	for (i = 0; i < INV_NUM; i++)
	{
		if (sigmaBaseNonZero[i])
		{
			for (k = 0; k < n; k++)
			{
				/////////////////////////////////
				if (TASK == 0)
				{
					sigma[k] = sigma[k] + a[i] * sigmaBase[idx];
					idx++;
				} 
				else if (TASK == 1)
				{
					if (mask[k] == 0)
					{
						sigma[k] = sigma[k] + a[i] * sigmaBase[idx];
						idx++;
					} 
					else
					{
						idx++;
					}
				}
				///////////////////////////////////
			}	
		}
		else
		{
			idx = idx + n;
		}
	}
	
	geneSigmaTV(O, sigmaTV, sigmaTVNonzero, p, q);
	if (sigmaTVNonzero)
	{
		for (k = 0; k < n; k++)
		{
			//sigma[k] += sigmaTV[k];
			sigma[k] = BETA_INV * sigma[k] + ALPHA_BV * sigmaTV[k]; 
		}
	}

	delete []inv;
	delete []sigmaBase;
	delete []sigmaTV;

	return;
}

void geneSigmaBase(DifImage_T &O, double *sigma, bool *sigmaBaseNonzero, int p, int q)
{
	int i, j, jStart;
	int n = O.getWidth() * O.getHeight();
	double *Ox = O.getVal(1, 0);
	double *Oy = O.getVal(0, 1);
	double *Oxx = O.getVal(2, 0);
	double *Oxy = O.getVal(1, 1);
	double *Oyy = O.getVal(0, 2);

	for (i = 0; i < INV_NUM; i++)
	{
		if (i == 0)
		{
			sigmaBaseNonzero[i] = false;
		}
		else
		{
			sigmaBaseNonzero[i] = true;
		}
	}

	// sigma_1:pq
	i = 1;
	jStart = i * n;
	if (p == 0 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 1;
		}
	}
	else
	{
		sigmaBaseNonzero[i] = false;
	}

	// sigma_2:pq
	i = 2;
	jStart = i * n;
	if (p == 1 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * Ox[j];
		}
	}
	else if (p == 0 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * Oy[j];
		}
	}
	else
	{
		sigmaBaseNonzero[i] = false;
	}

	// sigma_3:pq
	i = 3;
	jStart = i * n;
	if ((p == 2 && q == 0) || (p == 0 && q == 2))
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 1;
		}
	}
	else
	{
		sigmaBaseNonzero[i] = false;
	}

	// sigma_4:pq
	i = 4;
	jStart = i * n;
	if (p == 1 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * (Ox[j] * Oxx[j] + Oy[j] * Oxy[j]);
		}
	}
	else if (p == 0 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * (Ox[j] * Oxy[j] + Oy[j] * Oyy[j]);
		}
	}
	else if (p == 2 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = Ox[j] * Ox[j];
		}
	}
	else if (p == 1 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * Ox[j] * Oy[j];
		}
	}
	else if (p == 0 && q == 2)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = Oy[j] * Oy[j];
		}
	}
	else
	{
		sigmaBaseNonzero[i] = false;
	}

	// sigma_5:pq
	i = 5;
	jStart = i * n;
	if (p == 2 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * Oxx[j];
		}
	}
	else if (p == 1 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 4 * Oxy[j];
		}
	}
	else if (p == 0 && q == 2)
	{
		for (j = 0; j < n; j++)
		{
			sigma[jStart + j] = 2 * Oyy[j];
		}
	}
	else
	{
		sigmaBaseNonzero[i] = false;
	}

	return;
}

void geneSigmaTV(DifImage_T &O, double *sigmaTV, bool &sigmaTVNonzero, int p, int q)
{
	int j;
	int n = O.getWidth() * O.getHeight();
	double *Ox = O.getVal(1, 0);
	double *Oy = O.getVal(0, 1);
	double *Oxx = O.getVal(2, 0);
	double *Oxy = O.getVal(1, 1);
	double *Oyy = O.getVal(0, 2);
	
	double eps = 1e-6;
	double tmp1, tmp2, tmp3;
	sigmaTVNonzero = true;

	if (p == 0 && q == 0)
	{
		sigmaTVNonzero = false;
	}
	else if (p == 1 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] + Oy[j] * Oy[j];
			tmp2 = tmp1 * tmp1 * tmp1 + eps;
			tmp3 = sqrt(tmp2);
			sigmaTV[j] = 2 * (Ox[j] * Oyy[j] - Oy[j] * Oxy[j]) * tmp3 
				- 3 * Ox[j] * sqrt(tmp1) * (Ox[j] * Ox[j] * Oyy[j] 
				+ Oy[j] * Oy[j] * Oxx[j] - 2 * Ox[j] * Oy[j] * Oxy[j]);
			sigmaTV[j] = sigmaTV[j] / tmp2;
		}
	}
	else if ( p == 0 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] + Oy[j] * Oy[j];
			tmp2 = tmp1 * tmp1 * tmp1 + eps;
			tmp3 = sqrt(tmp2);
			sigmaTV[j] = 2 * (Oy[j] * Oxx[j] - Ox[j] * Oxy[j]) * tmp3 
				- 3 * Oy[j] * sqrt(tmp1) * (Ox[j] * Ox[j] * Oyy[j] 
				+ Oy[j] * Oy[j] * Oxx[j] - 2 * Ox[j] * Oy[j] * Oxy[j]);
			sigmaTV[j] = sigmaTV[j] / tmp2;
		}
	}
	else if ( p == 1 && q == 1)
	{
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] + Oy[j] * Oy[j];
			tmp2 = tmp1 * tmp1 * tmp1 + eps;
			tmp3 = sqrt(tmp2);
			sigmaTV[j] = (- 2 * Ox[j] * Oy[j] ) / tmp3;	
		}
	}
	else if (p == 2 && q == 0)
	{
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] + Oy[j] * Oy[j];
			tmp2 = tmp1 * tmp1 * tmp1 + eps;
			tmp3 = sqrt(tmp2);
			sigmaTV[j] = (Oy[j] * Oy[j]) / tmp3;
		}
	}
	else if (p == 0 && q == 2)
	{
		for (j = 0; j < n; j++)
		{
			tmp1 = Ox[j] * Ox[j] + Oy[j] * Oy[j];
			tmp2 = tmp1 * tmp1 * tmp1 + eps;
			tmp3 = sqrt(tmp2);
			sigmaTV[j] = (Ox[j] * Ox[j]) / tmp3;
		}
	}
		
	return;
}