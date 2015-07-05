
//************************************************************************
//*  Module Name: 
//*        implementation of image file input/output
//*  Abstract: 
//*        implements image read and write procedure
//*  Note: 
//*        none.
//* 
//*  Last Modified on: 2009-8-3 10:44:45 by Risheng Liu
//************************************************************************

#include "stdafx.h"
#include "GlobalApi.h"
#include "DifImage.h"
#include "VisCore.h"

//************************************************************************
//*  Function Name: imageRead
//*  Function Description: 
//*      read image from path.
//*  Arguments: 
//*      [OUT] : DifImage_T &inImage
//*      [IN] : string imageFile
//*      [IN] : int gapWidth
//* 
//*  Return Value: bool 
//*      bool value to judge if the read process is success.
//* 
//*  Last Modified on: 2009-8-3 10:40:14 by Risheng Liu
//************************************************************************

bool imageRead(DifImage_T &inImage, string imageFile, int gapWidth)
{
	bool success = true;
	int x, y;
	int width = inImage.getWidth();
	int height = inImage.getHeight();
	int n = width * height;
	unsigned char *pflt;
	double *val = new double[n];
	CVisByteImage image;
	if (!image.FReadFile(imageFile.data()))
	{
		success = false;
	}
	else
	{
		//copy pixel data to in
		for (y = 0; y < height; y++)
		{
			if (y >= gapWidth && y < height - gapWidth)
			{
				pflt = image.RowPointer(y - gapWidth);
				for (x = 0; x < width; x++)
				{
					if (x >= gapWidth && x < width - gapWidth)
					{
						val[x + y * width] = pflt[x - gapWidth] / (double)NORM_GRAY;
					}
					else
					{
						val[x + y * width] = 0;
					}
				}
			}
			else
			{
				for (x = 0; x < width; x++)
				{
					val[x + y * width] = 0;
				}
			}
		}
		inImage.setVal(val);
	}
	delete []val;

	return success;
}

//************************************************************************
//*  Function Name: imageWrite
//*  Function Description: 
//*      write image to path.
//*  Arguments: 
//*      [IN] : double *val
//*      [IN] : string imageFile
//*      [IN] : int width
//*      [IN] : int height
//*      [IN] : int gapWidth
//* 
//*  Return Value: bool 
//*      bool value to judge if the write process is success.
//* 
//*  Last Modified on: 2009-8-3 10:43:51 by Risheng Liu
//************************************************************************

bool imageWrite(double *val, string imageFile, int width, int height, int gapWidth)	//copy vector<double> val to pixel data
{
	bool success = true;
	int x, y;
	CVisByteImage image(width - 2 * gapWidth, height - 2 * gapWidth);
	unsigned char *pflt;
	double tmp;

	for (y = gapWidth; y < height - gapWidth; y++)
	{
		pflt = image.RowPointer(y - gapWidth);
		for (x = gapWidth; x < width - gapWidth; x++)
		{
			tmp = (int) (val[x + y * width] * NORM_GRAY + 0.5);
			if (tmp < 0)
			{
				pflt[x - 2] = 0;
			}
			else if (tmp > 255)
			{
				pflt[x - 2] = 255;
			}
			else
			{
				pflt[x - 2] = (unsigned char) tmp;
			}
		}
	}
	if (!image.FWriteFile(imageFile.data()))
	{
		success = false;
	}

	return success;
}