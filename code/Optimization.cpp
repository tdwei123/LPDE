
#include "StdAfx.h"
#include "GlobalApi.h"

extern ofstream logFile;

void minimize (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, GlobalParam_T &globalParam)
{
	int i;
	double phi, lastPhi, beta, gMod;

	int M = globalParam.getImageNum();
	int aLength = globalParam.getALength();
	double dt = globalParam.getDt();
	double eps = globalParam.getMinObjVal();
	int maxIter = globalParam.getMaxIter();
	int totalTime = int(TP / dt);

	double *g = new double[aLength];
	double *d = new double[aLength];
	double *da = new double[aLength];
	double *p = new double[aLength];
	double aTmp[INV_NUM];

	int counter = 0;

	//set initial values of a_ij(t)
	if (globalParam.getSelfInit())
	{
		init (inImage, maskImage, outImage, a, globalParam);
	}
	else
	{
		ifstream par;
		par.open("a.txt");
		if (par.good())
		{
			for (i = 0; i < aLength; i++)
			{
				par>>a[i];
			}
			par.close();
		}
		else
		{
			cout<<"Error of Reading a.txt!!"<<endl;
			return;
		}
	}

	//phi = eps + 1.0;
	phi = objVal(inImage, maskImage, outImage, a, a, 0, globalParam);
	cout<<"~"<<phi * NORM_GRAY * NORM_GRAY<<endl;
	logFile<<"~"<<phi * NORM_GRAY * NORM_GRAY<<endl;
	lastPhi = phi; //objVal(inImage, outImage, a, a, 0.0, globalParam);
	gMod = 1 / MIN_DOUBLE;
	while (lastPhi > eps && counter < maxIter)
	{
		//compute gradient
		cout<<"###grad"<<endl;
		grad(inImage, maskImage, outImage, a, g, globalParam);

		gMod = innerProduct(g, g, aLength);

		if (counter == 0)
		{
			for (i = 0; i < aLength; i++)
			{
				d[i] = -g[i];
			}
		}
		else
		{
			beta = 0.0;
			for (i = 0; i < aLength; i++)
			{
				beta = beta + g[i] * (g[i] + p[i]);
			}
			beta = beta / innerProduct(p, p, aLength);
			beta = (beta > 0) ? beta : 0;
			for (i = 0; i < aLength; i++)
			{
				d[i] = - g[i] + beta * d[i];
			}
			// restart
			/*if (innerProduct(d, g, aLength) >= 0)
			{
				cout<<"%Restart%"<<endl;
				for (i = 0; i < aLength; i++)
				{
					d[i] = - g[i];
				}
			}*/
		}
		for (i = 0; i < aLength; i++)
		{
			p[i] = - g[i];
		}

		// line minimization and update solution
		cout<<"Compute the linear minimization of the objective function!!"<<endl;

		phi = linMin(inImage, maskImage, outImage, a, d, da, globalParam);
		for (i = 0; i < aLength; i++)
		{
			a[i] = a[i] + da[i];
		}

		cout<<"~"<<phi * NORM_GRAY * NORM_GRAY<<endl;
		logFile<<"~"<<phi * NORM_GRAY * NORM_GRAY<<endl;

		for (i = 0; i < INV_NUM; i++)
		{
			aTmp[i] = a[i + aLength - INV_NUM];
		}
		restart(inImage, maskImage, outImage, a, globalParam);
		cout<<"restart:"<<endl;
		beta = phi;
		phi = objVal(inImage, maskImage, outImage, a, d, 0, globalParam);

		if (phi < beta)
		{
			setZeroValue(d, aLength);
		}
		else
		{
			for (i = 0; i < INV_NUM; i++)
			{
				a[i + aLength - INV_NUM] = aTmp[i];
			}
			phi = beta;
		}
		if (fabs(lastPhi - phi) < 1e-7 * phi)
		{
			break;
		}
		lastPhi = phi;

		counter++;
		cout<<"Counter: "<<counter<<endl;

		/*if (counter % 10 == 0)
		{
			cout<<"Counter: "<<counter<<endl;
		}*/
	}

	delete []d;
	delete []da;
	delete []g;
	delete []p;

	return;
}

void grad (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *g, GlobalParam_T &globalParam)
{
	int i, k, s, t, timeInd, idx, ind;
	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;

	int M = globalParam.getImageNum();
	int aLength = globalParam.getALength();
	double dt = globalParam.getDt();
	double *pLamta = globalParam.getLamta();
	int totalTime = int(TP / dt);

	double *mask;

	double *pLagMul = new double[totalTime * n];
	DifImage_T *pGen = new DifImage_T[(totalTime + 1)];
	for (i = 0; i <= totalTime; i++)
	{
		pGen[i].setArea(width, height);
	}
	double *inv = new double[INV_NUM * n];

	for (i = 0; i < aLength; i++)
	{
		g[i] = 0;
	}

	for (s = 0; s < M; s++)
	{
		mask = maskImage[s].getVal(0, 0);

		// get lagrange multiplier phi_m(t)
		pdeAdjSolver(inImage[s], maskImage[s], outImage[s], pGen, pLagMul, a, dt);
		timeInd = 0;
		ind = 0;
		for (t = 0; t < totalTime; t++)
		{
			geneInv(pGen[t], inImage[s], inv);
			idx = 0;
			for (i = 0; i < INV_NUM; i++)
			{
				for (k = 0; k < n; k++)
				{
					/////////////////////////////////////
					if (TASK == 0)
					{
						g[ind] = g[ind] - pLagMul[timeInd + k] * inv[idx];
						idx++;
					}
					else if (TASK == 1)
					{
						if (mask[k] == 0)
						{
							g[ind] = g[ind] - pLagMul[timeInd + k] * inv[idx];
							idx++;
						}
						else
						{
							idx++;
						}
					}
					////////////////////////////////////
				}
				ind++;
			}
			timeInd = timeInd + n;
		}
	}

	for (i = 0; i < aLength; i++)
	{
		g[i] = (BETA_INV * g[i] / n + pLamta[i % INV_NUM] * a[i]) * dt;
	}

	delete []pLagMul;
	delete []pGen;
	delete []inv;
	
	return;
}

double linMin (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *d, double *da, GlobalParam_T &globalParam)
{
	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;
	int i;
	double w[4], wn, lim, p[4], ene, weps;

	int aLength = globalParam.getALength();

	/*weps = 0.0;
	for (i = 0; i < aLength; i++)
	{
		if (fabs(d[i]) > weps)
		{
			weps = fabs(d[i]);
		}
	}*/
	weps = 1.0e-7 / (innerProduct(d, d, aLength) + MIN_DOUBLE);

	//bracket
	cout<<"$$$bracket"<<endl;

	w[0] = 0.0;
	p[0] = objVal(inImage, maskImage, outImage, a, d, w[0], globalParam);
	w[1] = weps;   // intial bracket for bracket search
	p[1] = objVal(inImage, maskImage, outImage, a, d, w[1], globalParam);
	if (p[0] < p[1])
	{
		wn = w[0];
		w[0] = w[1];
		w[1] = wn;
	}
	w[2] = w[1] + 1.618 * (w[1] - w[0]);
	p[2] = objVal(inImage, maskImage, outImage, a, d, w[2], globalParam);
	while (p[1] > p[2] && w[1] != w[2])
	{
		w[3] = 0.5 * ((w[1] * w[1] - w[2] * w[2]) * p[0] + (w[2] * w[2] - w[0] * w[0]) * p[1]
			+ (w[0] * w[0] - w[1] * w[1]) * p[2]) /
			((w[1] - w[2]) * p[0] + (w[2] - w[0]) * p[1] + (w[0] - w[1]) * p[2] + MIN_DOUBLE);
		lim = w[1] + 1000 * (w[2] - w[1]);
		cout<<"-----------exterpolation: "<<w[0]<<" "<<w[1]<<" "<<w[2]<<" "<<w[3]<<endl;
		if ((w[3] - w[2]) * (w[1] - w[3]) > 0)    //w[3] is between w[1] and w[2]
		{
			p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
			cout<<"-----------exterpolation: "<<p[2]<<" "<<p[3]<<"w: "<<w[2]<<" "<<w[3]<<endl;
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
			p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
		}
		else if ((w[2] - w[3]) * (w[3] - lim) > 0) //w[3] is between w[2] and lim
		{
			p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
			if (p[3] < p[2])
			{
				w[1] = w[2];
				p[1] = p[2];
				w[2] = w[3];
				p[2] = p[3];
				w[3] = w[2] + 1.618 * (w[2] - w[1]);
				p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
			}
		}
		else if ((w[3] - lim) * (lim - w[2]) >= 0)
		{
			w[3] = lim;
			p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
		}
		else
		{
			w[3] = w[2] + 1.618 * (w[2] - w[1]);
			p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
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
	cout<<"@@@search"<<endl;
	w[3] = w[2];              //bracket (w[0], w[3])
	w[1] = w[3] - 0.618 * (w[3] - w[0]);
	w[2] = w[0] + 0.618 * (w[3] - w[0]);
	p[1] = objVal(inImage, maskImage, outImage, a, d, w[1], globalParam);
	p[2] = objVal(inImage, maskImage, outImage, a, d, w[2], globalParam);
	while(fabs(w[3] - w[0]) > 0.01 * (weps + fabs(w[1]) + fabs(w[2]))) //0.001 * (fabs(w[1]) + fabs(w[2]) + 0.1))
	{
		cout<<"***********golden section: "<<w[0]<<" "<<w[3]<<endl;
		//cout<<p[1]<<" "<<p[2]<<endl;
		if (p[1] > p[2])
		{
			w[0] = w[1];
			w[1] = w[2];
			w[2] = w[0] + 0.618 * (w[3] - w[0]);
			p[1] = p[2];
			p[2] = objVal(inImage, maskImage, outImage, a, d, w[2], globalParam);
		}
		else //if (p[1] < p[2])
		{
			w[3] = w[2];
			w[2] = w[1];
			w[1] = w[3] - 0.618 * (w[3] - w[0]);
			p[2] = p[1];
			p[1] = objVal(inImage, maskImage, outImage, a, d, w[1], globalParam);
		}
	}
	p[0] = objVal(inImage, maskImage, outImage, a, d, w[0], globalParam);
	p[3] = objVal(inImage, maskImage, outImage, a, d, w[3], globalParam);
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

	cout<<"linmin coefficient: "<<wn<<endl;
	//update solution
	for (i = 0; i < aLength; i++)
	{
		da[i] = wn * d[i];
	}

	/*if (fabs(wn) < 1e-6 * weps)
	{
		return -1;
	}*/
	return ene;
}

double objVal (DifImage_T *inImage, DifImage_T *maskImage, DifImage_T *outImage, double *a, double *d, double w, GlobalParam_T &globalParam)
{
	int i, k;
	double phi;

	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;
	int M = globalParam.getImageNum();
	int aLength = globalParam.getALength();
	double dt = globalParam.getDt();
	double *pLamta = globalParam.getLamta();
	int totalTime = int(TP / dt);

	double *pGeneImageData1 = new double[n];
	double tmp;

	double *aTmp = new double[aLength];
	for (i = 0; i < aLength; i++)
	{
		aTmp[i] = a[i] + w * d[i];
	}

	phi = 0.0;
	for (i = 0; i < M; i++)
	{
		// solve and do not need to store O because we only need O(1) and dO/dt = L(I, O, t), store phi
		pdeSolver(inImage[i], maskImage[i], pGeneImageData1, aTmp, dt);
		for (k = 0; k < n; k++)
		{
			/*if(fabs(outImage[i].getVal(0, 0, k)) > 10 || fabs(pGeneImageData1[k]) > 10)
			{
				cout<<outImage[i].getVal(0, 0, k);
			}*/
			tmp = outImage[i].getVal(0, 0, k) - pGeneImageData1[k];
			phi = phi + tmp * tmp;
		}
	}
	phi = phi / n;

	tmp = 0.0;
	for (i = 0; i < aLength; i++)
	{
		tmp = tmp + pLamta[i % INV_NUM] * aTmp[i] * aTmp[i];
	}
	tmp = tmp * dt;
	phi = 0.5 * (phi + tmp);

	delete []aTmp;
	delete []pGeneImageData1;

	//logFile<<"objective function value: "<<phi * NORM_GRAY * NORM_GRAY<<endl;
	cout<<"objective function value: "<<phi * NORM_GRAY * NORM_GRAY<<endl;

	return phi;
}
