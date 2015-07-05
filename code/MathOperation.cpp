
//************************************************************************
//*  Module Name: 
//*        mathematical operators and optimal algorithms for LPDE
//*  Abstract: 
//*        
//*  Note: 
//*        none.
//* 
//*  Last Modified on: 2009-8-3 10:59:45 by Risheng Liu
//************************************************************************

#include "StdAfx.h"
#include "GlobalApi.h"

//************************************************************************
//*  Function Name: multiply
//*  Function Description: 
//*      matrix multiply vector
//*  Arguments: 
//*      [IN] : double *A - input matrix
//*      [IN] : double *b - input vector
//*      [OUT] : double *Ab - output vector
//*      [OUT] : int n - length of vector
//* 
//*  Return Value: void 
//*      none.
//* 
//*  Last Modified on: 2009-8-3 11:35:49 by Risheng Liu
//************************************************************************

void multiply (double *A, double *b, double *Ab, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		Ab[i] = 0.0;
		for (j = 0; j < n; j++)
		{
			Ab[i] = Ab[i] + A[i * n + j] * b[j];
		}
	}
	return;
}

//************************************************************************
//*  Function Name: innerProduct
//*  Function Description: 
//*      vector multiply vector
//*  Arguments: 
//*      [IN] : double *a
//*      [IN] : double *b
//*      [IN] : int n
//* 
//*  Return Value: double 
//*      the inner product of two vector.
//* 
//*  Last Modified on: 2009-8-3 12:40:20 by Risheng Liu
//************************************************************************

double innerProduct (double *a, double *b, int n)
{
	int i;
	double sum;
	sum = 0.0;
	for (i = 0; i < n; i++)
	{
		sum = sum + a[i] * b[i];
	}
	return sum;
}