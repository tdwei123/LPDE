// LPdeO.cpp : Defines the entry point for the console application.
//


#include "StdAfx.h"
#include "windows.h"

#include "GlobalApi.h"

ofstream logFile;

int main(int argc, char* argv[])
{
	GlobalParam_T globalParam;
	logFile.open("_log.txt", ios_base::app);
	int M;
	double dt, eps;
	int width, height, gapWidth;
	double lamta[INV_NUM];
	int maxIter;
	bool selfInit = true;
	int i;

	time_t t;
	t = time(NULL);
	logFile<<"========================================"<<endl<<ctime(&t)<<endl;

	ifstream dirFile, paramFile;
	if (TASK == 0)
	{
		dirFile.open("dir_denoise.txt");
		paramFile.open("param_denoise.txt");
	}
	else if (TASK == 1)
	{
		dirFile.open("dir_inpaint.txt");
		paramFile.open("param_inpaint.txt");
	}
	
	
	if (!dirFile.good()||!paramFile.good())	
	{
		cout<<"~Can not open the files!!"<<endl
			<<"~Please Enter 0 to End The Program!!"<<endl;
	    while (getchar() != '0')
		{
		}
		return 0;
	}
	
	paramFile>>lamta[0];
	for (i = 1; i < INV_NUM; i++)
	{
		lamta[i] = lamta[0];
	}
	paramFile>>dt;
	paramFile>>M;
	paramFile>>eps;
	paramFile>>width;
	paramFile>>height;
	paramFile>>gapWidth;
	paramFile>>maxIter;
	paramFile>>selfInit;

	globalParam.setLamta(lamta);
	globalParam.setDt(dt);
	globalParam.setImageNum(M);
	globalParam.setGapWidth(gapWidth);
	globalParam.setMinObjVal(eps * M / NORM_GRAY / NORM_GRAY);
	globalParam.setHeight(height + 2 * gapWidth);
	globalParam.setWidth(width + 2 * gapWidth);
	globalParam.setMaxIter(maxIter);
	globalParam.setSelfInit(selfInit);

	char c[100];

	string rootDir, inFilePath, outFilePath, geneFilePath, testFilePath,
		taskType[4], inMaskPath, testMaskPath, mask[2];

	WIN32_FIND_DATA findData;
	HANDLE hdlFind;

	dirFile.getline(c, 50);
	rootDir = c;
	dirFile.getline(c, 50);
	taskType[0] = c;
	dirFile.getline(c, 50);
	taskType[1] = c;
	dirFile.getline(c, 50);
	taskType[2] = c;
	dirFile.getline(c, 50);
	taskType[3] = c;

	dirFile.getline(c, 50);
	mask[0] = c;
	dirFile.getline(c, 50);
	mask[1] = c;

	inFilePath = rootDir + "\\" + taskType[0];
	outFilePath = rootDir + "\\" + taskType[1];
	geneFilePath = rootDir + "\\" + taskType[3];

	hdlFind = ::FindFirstFile(geneFilePath.data(), &findData);
	if (hdlFind == INVALID_HANDLE_VALUE)
	{
		int k = geneFilePath.length();
	    for (i = 0; i < k; i++)
		{
			if (geneFilePath[i] == '\\')
			{
				hdlFind = ::FindFirstFile(geneFilePath.substr(0, i).data(), &findData);
				if (hdlFind == INVALID_HANDLE_VALUE)
				{
					CreateDirectory(geneFilePath.substr(0, i + 1).data(), NULL);
				}
			}
		}
		CreateDirectory(geneFilePath.data(), NULL);
	}
	testFilePath = rootDir + "\\" + taskType[2];

	if (TASK == 1)
	{
		inMaskPath = rootDir + "\\" + mask[0];
		testMaskPath = rootDir + "\\" + mask[1];
	}
	else if (TASK == 0)
	{
		inMaskPath = inFilePath;
		testMaskPath = testFilePath;
	}

	int aLength = globalParam.getALength();
	double *a = new double[aLength];

	//initA(a, aLength);
	train(inFilePath, outFilePath, geneFilePath, inMaskPath, taskType[3], a, globalParam);
	test(testFilePath, geneFilePath, testMaskPath, a, globalParam);

	delete []a;

	dirFile.close();
	paramFile.close();
	logFile.close();

	cout<<"Please Enter 1 to End The Program!"<<endl;
	while (getchar() != '1')
	{
	}
		
	return 0;
}


//************************************************************************
//*  Function Name: train 
//*  Function Description: 
//*      learn coefficient function from image pairs
//*  Arguments: 
//*      [IN] : string inFilePath
//*      [IN] : string outFilePath
//*      [IN] : string geneFilePath
//*      [IN] : string taskType
//*      [OUT] : double *a
//*      [IN] : GlobalParam_T &globalParam
//* 
//*  Return Value: void 
//*      none.
//* 
//*  Last Modified on: 2009-8-3 10:47:48 by Risheng Liu
//************************************************************************

void train (string inFilePath, string outFilePath, string geneFilePath, string inMaskPath, string taskType, double *a, GlobalParam_T &globalParam)
{
	BOOL bWorking;
	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;
	int M = globalParam.getImageNum();
	int aLength = globalParam.getALength();
	double dt = globalParam.getDt();
	int gapWidth = globalParam.getGapWidth();
	DifImage_T *inImage = new DifImage_T[M];
	DifImage_T *outImage = new DifImage_T[M];

	DifImage_T *maskImage = new DifImage_T[M];
	int i;

	for (i = 0; i < M; i++)
	{
		inImage[i].setArea(width, height);
		maskImage[i].setArea(width, height);
		outImage[i].setArea(width, height);
	}

	string inFileName;
	string outFileName;
	string maskFileName;

	string fileName;
	string aFName;

	WIN32_FIND_DATA findData;
	fileName = inFilePath + "\\*.bmp";
	HANDLE hdlFind = ::FindFirstFile(fileName.data(), &findData);
	for (i = 0; i < M && bWorking; i++)
	{
		do 
		{
			if (!bWorking)
			{
				break;
			}
			fileName = findData.cFileName;

			inFileName = inFilePath + "\\" + fileName;
			maskFileName = inMaskPath + "\\" + fileName; 
			outFileName = outFilePath + "\\" + fileName;
			cout<<"In File Name:"<<inFileName<<endl;
			if (TASK)
			{
				cout<<"Mask File Name:"<<maskFileName<<endl;
			}
			cout<<"Out File Name:"<<outFileName<<endl;
			imageRead (inImage[i], inFileName, gapWidth);
			imageRead(maskImage[i], maskFileName, gapWidth);
			bWorking = ::FindNextFile (hdlFind, &findData);

		} while(!imageRead (outImage[i], outFileName, gapWidth));			
	}
	::FindClose(hdlFind);

	if (i < M)
	{
		M = i;
		globalParam.setImageNum(M);
	}

	time_t t;
	t = time(NULL);
	minimize (inImage, maskImage, outImage, a, globalParam);
	logFile<<"Time Cost:"<<time(NULL) - t<<endl;
	cout<<"Time Cost:"<<time(NULL) - t<<endl;

	ofstream aFile;
	aFName = geneFilePath + "\\a.txt";
	aFile.open(aFName.data(), ios_base::app);
	aFile<<endl<<"============================================="<<endl<<ctime(&t)<<endl;
	aFile<<"TaskType: "<<taskType<<endl;

	aFile<<"a=: ";
	for (i = 0; i < aLength; i++)
	{
		if (i % INV_NUM == 0)
		{
			aFile<<endl;
		}
		aFile<<a[i]<<" ";
	}
	aFile.close();

	delete []inImage;
	delete []outImage;
	return;
}

//************************************************************************
//*  Function Name: test 
//*  Function Description: 
//*      test the learnt PDE and save the result image
//*  Arguments: 
//*      [IN] : string testFilePath
//*      [IN] : string geneFilePath
//*      [IN] : double *a
//*      [IN] : GlobalParam_T &globalParam
//* 
//*  Return Value: void 
//*      use learnt PDE to solve the version task
//* 
//*  Last Modified on: 2009-8-3 10:51:20 by Risheng Liu
//************************************************************************

void test (string testFilePath, string geneFilePath, string testMaskPath, double *a, GlobalParam_T &globalParam)
{
	string fileName, testFileName, geneFileName, maskFileName;
	int width = globalParam.getWidth();
	int height = globalParam.getHeight();
	int n = width * height;
	double dt = globalParam.getDt();
	int gapWidth = globalParam.getGapWidth();
	DifImage_T testImage;
	DifImage_T maskImage;
	testImage.setArea(width, height);
	maskImage.setArea(width, height);// important
	double *pGeneImageData1 = new double[n];
	
	testFileName = testFilePath + "\\*.bmp";
	maskFileName = testMaskPath + "\\*.bmp";
	WIN32_FIND_DATA findData;
	HANDLE hdlFind = ::FindFirstFile(testFileName.data(), &findData);
	BOOL bWorking = (hdlFind != INVALID_HANDLE_VALUE);
	while (bWorking)
	{
		fileName = findData.cFileName;
		testFileName = testFilePath + "\\" + fileName;
		maskFileName = testMaskPath + "\\" + fileName;
		geneFileName = geneFilePath + "\\" + fileName;
		cout<<"test file name: "<<testFileName<<endl;
		imageRead(testImage, testFileName, gapWidth);
		imageRead(maskImage, maskFileName, gapWidth);
		pdeSolver(testImage, maskImage, pGeneImageData1, a, dt);
		imageWrite(pGeneImageData1, geneFileName, width, height, gapWidth);

		// test
		//double *tmp = testImage.getVal(0,0);
		//GaussFilter(tmp, width, height, 7, 5);
		//imageWrite(tmp, geneFileName, width, height, gapWidth);
		// test

		bWorking = ::FindNextFile(hdlFind, &findData);
	}
	::FindClose(hdlFind);
	delete []pGeneImageData1;
}


void showStepResult(double *val)
{
	///////Test The Step Result/////////////
	int gapWidth = 2;
	string geneFileName =  "test.bmp";
	double width = 485;
	double height = 325;
	imageWrite(val, geneFileName, width, height, gapWidth);
	////////////////////////////////////////
}

void initA(double *a, double aLength)
{
	int i;
	ifstream par;
	par.open("a.txt");
	if (par.good())
	{
		for (i = 0; i < aLength; i++)
		{
			par>>a[i];
		}
		par.close();
		cout<<"Initial from a.txt!"<<endl;
	}
	else
	{
		cout<<"Error of Reading a.txt!!"<<endl;
		return;
	}

}
