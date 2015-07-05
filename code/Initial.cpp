
#include "StdAfx.h"
#include "GlobalApi.h"

void init (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam)
{
	// O_m(k dt) -> a(k dt) -> O_m((k+1)dt)
	int s, t;
	double *pA;

	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int M = globalParam.getImageNum();
	int aLength = globalParam.getALength();
	double dt = globalParam.getDt();
	double eps = globalParam.getMinObjVal();
	int totalTime = int(TP / dt);
	int gapWidth = globalParam.getGapWidth();

	DifImage_T *pGen = new DifImage_T[M];
	DifImage_T nextGen;
	nextGen.setArea(width, height);

	for (s = 0; s < M; s++)
	{
		pGen[s].setArea(width, height);
		pGen[s].setVal(inImage[s].getVal(0, 0));
	}

	setZeroValue(a, aLength); // assign a_i to zero
	pA = a;
	for (t = 0; t < totalTime; t++)
	{
		stepInit(pGen, inImage, maskImage, outImage, pA, t, globalParam);
		if (t < totalTime - 1)
		{
			for (s = 0; s < M; s++)
			{
				pdeStep(pGen[s], nextGen, inImage[s], maskImage[s], pA, dt);
				pGen[s].setVal(nextGen.getVal(0, 0));
			}
		}
		pA = pA + INV_NUM;
	}

	/////////////////////////////
	cout<<"a: ";
	for (t = 0; t < INV_NUM * totalTime; t++)
	{
		if (t % INV_NUM == 0)
		{
			cout<<endl;
		}
		cout<<a[t]<<" ";
	}
	cout<<endl;
	////////////////////////////

	delete []pGen;
}

void restart (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam)
{
	// O_m(k dt) -> a(k dt) -> O_m((k+1)dt)
	int s, t;
	double *pA;

	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int M = globalParam.getImageNum();
	double dt = globalParam.getDt();
	int totalTime = int(TP / dt);

	DifImage_T *pGen = new DifImage_T[M];
	DifImage_T nextGen;
	nextGen.setArea(width, height);

	for (s = 0; s < M; s++)
	{
		pGen[s].setArea(width, height);
		pGen[s].setVal(inImage[s]);
	}
	pA = a;
	for (t = 0; t < totalTime; t++)
	{
		if (t == totalTime - 1)
		{
			stepInit(pGen, inImage, maskImage, outImage, pA, t, globalParam);
		}
		else
		{
			for (s = 0; s < M; s++)
			{
				pdeStep(pGen[s], nextGen, inImage[s], maskImage[s], pA, dt);
				pGen[s].setVal(nextGen);
			}
		}
		pA = pA + INV_NUM;
	}

	delete []pGen;
}

void stepInit (DifImage_T *pGen, DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, int t, GlobalParam_T &globalParam)
{
	int i, j, k, s, iStart, jStart;
	double c = 0.0;

	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;
	int M = globalParam.getImageNum();
	double dt = globalParam.getDt();

	double *inv = new double[INV_NUM * n];
	double *bv = new double[2 * n];
	
	double *dO = new double[n];
	double *O, *Out;
	double G[INV_NUM * INV_NUM], g[INV_NUM], tmp;

	double *mask;

	for (i = 0; i < INV_NUM; i++)
	{
		g[i] = 0.0;
		for (j = 0; j < INV_NUM; j++)
		{
			G[i * INV_NUM + j] = 0.0;
		}
	}
	for (s = 0; s < M; s++)
	{
		O = pGen[s].getVal(0, 0);
		Out = outImage[s].getVal(0, 0);

		mask = maskImage[s].getVal(0, 0);

		for (k = 0; k < n; k++)
		{
			dO[k] = (Out[k] - O[k]) / (TP - t * dt);
			c = c + dO[k] * dO[k];
		}

		geneInv(pGen[s], inImage[s], inv);
		geneBV(pGen[s], inImage[s], maskImage[s], bv);
		
		iStart = 0;
		for (i = 0; i < INV_NUM; i++)
		{
			jStart = iStart;
			for (j = i; j < INV_NUM; j++)
			{
				tmp = 0.0;
				for (k = 0; k < n; k++)
				{
					//tmp = tmp + inv[iStart + k] * inv[jStart + k]; // normalization later
					//////////////////////////////////////////
					if (TASK == 0)
					{
						//if (k % width > 1 && k % width < width - 2)
						tmp = tmp + inv[iStart + k] * inv[jStart + k]; // normalization later
					}
					else if (TASK == 1)
					{
						if (mask[k] == 0)
						{
							//if (k % width > 1 && k % width < width - 2)
							tmp = tmp + inv[iStart + k] * inv[jStart + k]; // normalization later
						} 
					}
					//////////////////////////////////////////
				}
				//G[i * INV_NUM + j] = G[i * INV_NUM + j] + tmp;
				G[i * INV_NUM + j] = G[i * INV_NUM + j] + BETA_INV * tmp;
				jStart = jStart + n;
			}
			tmp = 0.0;
			for (k = 0; k < n; k++)
			{
				//if (k % width > 1 && k % width < width - 2)
				tmp = tmp + inv[iStart + k] * (dO[k] - ALPHA_BV * (bv[k] + bv[n + k]));// bv[k]--tv, bv[k + n]--sgn
				//tmp = tmp + inv[iStart + k] * dO[k];
			}
			g[i] = g[i] + BETA_INV * tmp;
			iStart = iStart + n;
		}
	}
	for (i = 0; i < INV_NUM; i++)
	{
		for (j = 0; j < i; j++)
		{
			G[i * INV_NUM + j] = G[j * INV_NUM + i];  // G is symmetric
		}
		g[i] = g[i] / n;  // normalize g
	}
	for (i = 0; i < INV_NUM * INV_NUM; i++)
	{
		G[i] = G[i] / n;  // normalize G
	}
	c = c / n;

	qcqp(G, g, c, a);

	delete []inv;
	delete []bv;
	delete []dO;

	return;
}

void qcqp (double *F, double *b, double c, double *a)
{
	int i;

	for (i = 0; i < INV_NUM; i++)
	{
		a[i] = 0;
	}

	double g[INV_NUM], p[INV_NUM], d[INV_NUM];
	double r = 1e-6;
	int counter = 0;
	double beta, phi;
	double eps = 1e-6;
	double lastPhi = eps + 1.0;
	while (lastPhi > eps)
	{
		//compute gradient
		//cout<<"###grad"<<endl;
		multiply(F, a, g, INV_NUM);
		for (i = 0; i < INV_NUM; i++)
		{
			g[i] = g[i] - b[i];
		}


		if (counter == 0)
		{
			for (i = 0; i < INV_NUM; i++)
			{
				d[i] = -g[i];
			}
		}
		else
		{
			beta = 0.0;
			for (i = 0; i < INV_NUM; i++)
			{
				beta = beta + g[i] * (g[i] + p[i]);
			}
			beta = beta / innerProduct(p, p, INV_NUM);
			beta = (beta > 0) ? beta : 0;
			for (i = 0; i < INV_NUM; i++)
			{
				d[i] = - g[i] + beta * d[i];
			}
		}
		for (i = 0; i < INV_NUM; i++)
		{
			p[i] = - g[i];
		}

		//line minimization and update solution
		phi = lineFunW(F, b, c, a, d);

		if (fabs(lastPhi - phi) < 1e-6 * phi)
		{
			break;
		}
		lastPhi = phi;

		//cout<<phi<<endl;

		counter++;
	}

	cout<<"~"<<phi<<endl;

	return;
}

double lineFunW (double *F, double *b, double c, double *a, double *d)
{
	int i;
	double w[4], wn, lim, p[4], ene, weps;

	/*weps = 0.0;
	for (i = 0; i < aLength; i++)
	{
		if (fabs(d[i]) > weps)
		{
			weps = fabs(d[i]);
		}
	}*/
	weps = 1.0e-7 / (innerProduct(d, d, INV_NUM) + MIN_DOUBLE);

	//bracket

	w[0] = 0.0;
	p[0] = funW (F, b, c, a, d, w[0]);
	w[1] = weps;   // intial bracket for bracket search
	p[1] = funW (F, b, c, a, d, w[1]);
	if (p[0] < p[1])
	{
		wn = w[0];
		w[0] = w[1];
		w[1] = wn;
	}
	w[2] = w[1] + 1.618 * (w[1] - w[0]);
	p[2] = funW (F, b, c, a, d, w[2]);
	while (p[1] > p[2] && w[1] != w[2] && fabs(w[2]) < 1e5)
	{
		w[3] = 0.5 * ((w[1] * w[1] - w[2] * w[2]) * p[0] + (w[2] * w[2] - w[0] * w[0]) * p[1]
			+ (w[0] * w[0] - w[1] * w[1]) * p[2]) /
			((w[1] - w[2]) * p[0] + (w[2] - w[0]) * p[1] + (w[0] - w[1]) * p[2] + MIN_DOUBLE);
		lim = w[1] + 1000 * (w[2] - w[1]);
		//cout<<"-----------exterpolation: "<<w[0]<<" "<<w[1]<<" "<<w[2]<<" "<<w[3]<<endl;
		if ((w[3] - w[2]) * (w[1] - w[3]) > 0)    //w[3] is between w[1] and w[2]
		{
			p[3] = funW (F, b, c, a, d, w[3]);
			//cout<<"-----------exterpolation: "<<p[2]<<" "<<p[3]<<"w: "<<w[2]<<w[3]<<endl;
			if (p[3] < p[2])                      //get a minimum between w[1] and w[2]
			{
				w[0] = w[1];
				break;
			}
			else if (p[3] > p[1])                 //get a minimum between w[0] and w[3]
			{
				w[2] = w[3];
				break;
			}
			w[3] = w[2] + 1.618 * (w[2] - w[1]);
			p[3] = funW (F, b, c, a, d, w[3]);
		}
		else if ((w[2] - w[3]) * (w[3] - lim) > 0) //w[3] is between w[2] and lim
		{
			p[3] = funW (F, b, c, a, d, w[3]);
			if (p[3] < p[2])
			{
				w[1] = w[2];
				p[1] = p[2];
				w[2] = w[3];
				p[2] = p[3];
				w[3] = w[2] + 1.618 * (w[2] - w[1]);
				p[3] = funW (F, b, c, a, d, w[3]);
			}
		}
		else if ((w[3] - lim) * (lim - w[2]) >= 0)
		{
			w[3] = lim;
			p[3] = funW (F, b, c, a, d, w[3]);
		}
		else
		{
			w[3] = w[2] + 1.618 * (w[2] - w[1]);
			p[3] = funW (F, b, c, a, d, w[3]);
		}
		w[0] = w[1];
		p[0] = p[1];
		w[1] = w[2];
		p[1] = p[2];
		w[2] = w[3];
		p[2] = p[3];
	}
	//output bracket (w[0], w[2])

	//search
	//cout<<"@@@search"<<endl;
	w[3] = w[2];              //bracket (w[0], w[3])
	w[1] = w[3] - 0.618 * (w[3] - w[0]);
	w[2] = w[0] + 0.618 * (w[3] - w[0]);
	p[1] = funW (F, b, c, a, d, w[1]);
	p[2] = funW (F, b, c, a, d, w[2]);
	while(fabs(w[3] - w[0]) > 0.01 * (weps + fabs(w[1]) + fabs(w[2]))) //0.001 * (fabs(w[1]) + fabs(w[2]) + 0.1))
	{
		//cout<<"***********golden section: "<<w[0]<<" "<<w[3]<<endl;
		//cout<<p[1]<<" "<<p[2]<<endl;
		if (p[1] > p[2])
		{
			w[0] = w[1];
			w[1] = w[2];
			w[2] = w[0] + 0.618 * (w[3] - w[0]);
			p[1] = p[2];
			p[2] = funW (F, b, c, a, d, w[2]);
		}
		else //if (p[1] < p[2])
		{
			w[3] = w[2];
			w[2] = w[1];
			w[1] = w[3] - 0.618 * (w[3] - w[0]);
			p[2] = p[1];
			p[1] = funW (F, b, c, a, d, w[1]);
		}
	}
	p[0] = funW (F, b, c, a, d, w[0]);
	p[3] = funW (F, b, c, a, d, w[3]);
	wn = w[0];
	ene = p[0];
	for (i = 1; i <= 3; i++)
	{
		if (p[i] < ene)
		{
			wn = w[i];
			ene = p[i];
		}
	}

	for (i = 0; i < INV_NUM; i++)
	{
		a[i] = a[i] + wn * d[i];
	}

	return ene;
}

double funW (double *F, double *b, double c, double *a, double *d, double w)
{
	int i;
	double phi;
	double r = 1e-6;

	double aTmp[INV_NUM], p[INV_NUM];
	for (i = 0; i < INV_NUM; i++)
	{
		aTmp[i] = a[i] + w * d[i];
	}

	multiply(F, aTmp, p, INV_NUM);

	phi = (innerProduct(aTmp, p, INV_NUM) + c) / 2 - innerProduct(aTmp, b, INV_NUM);

	return phi;
}