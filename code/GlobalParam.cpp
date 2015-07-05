
#include "StdAfx.h"
#include "GlobalApi.h"
#include "GlobalParam.h"

GlobalParam_T::GlobalParam_T (void)
{
	for (int i= 0; i < INV_NUM; i++)
	{
		i_pLamta[i] = 0;
	}
	i_dt = 0.0;
	i_aLength = 0;
	i_imageNum = 0;
	i_minObjVal = 0.0;
	i_gapWidth = 0;
	i_width = 0;
	i_height = 0;
	maxIter = 0;
	selfInit = 1;
}

GlobalParam_T::GlobalParam_T (GlobalParam_T &p)
{
	for (int i = 0; i < INV_NUM; i++)
	{
		i_pLamta[i] = p.i_pLamta[i];
	}
	i_dt = p.i_dt;
	i_aLength = p.i_aLength;
	i_imageNum = p.i_imageNum;
	i_minObjVal = p.i_minObjVal;
	i_gapWidth = p.i_gapWidth;
	i_width = p.i_width;
	i_height = p.i_height;
	maxIter = p.maxIter;
	selfInit = p.selfInit;
}

GlobalParam_T::~GlobalParam_T (void)
{
}

void GlobalParam_T::setLamta (double *pLamta)
{
	for (int i = 0; i < INV_NUM; i++)
	{
		i_pLamta[i] = pLamta[i];
	}
}

void GlobalParam_T::setDt (double dt)
{
	i_dt = dt;
	i_aLength = INV_NUM * int(TP / dt);
}

void GlobalParam_T::setImageNum (int imageNum)
{
	i_imageNum = imageNum;
}

void GlobalParam_T::setMinObjVal (double eps)
{
	i_minObjVal = eps;
}

void GlobalParam_T::setGapWidth (int gapWidth)
{
	i_gapWidth = gapWidth;
}

void GlobalParam_T::setWidth (int width)
{
	i_width = width;
}

void GlobalParam_T::setHeight (int height)
{
	i_height = height;
}

double *GlobalParam_T::getLamta (void)
{
	return i_pLamta;
}

double GlobalParam_T::getDt (void)
{
	return i_dt;
}

int GlobalParam_T::getALength (void)
{
	return i_aLength;
}

int GlobalParam_T::getImageNum (void)
{
	return i_imageNum;
}

double GlobalParam_T::getMinObjVal (void)
{
	return i_minObjVal;
}

int GlobalParam_T::getGapWidth (void)
{
	return i_gapWidth;
}

int GlobalParam_T::getWidth (void)
{
	return i_width;
}

int GlobalParam_T::getHeight (void)
{
	return i_height;
}
