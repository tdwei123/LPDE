
#pragma once

#define  INV_NUM 6

class GlobalParam_T {
public:
	GlobalParam_T (void);
	GlobalParam_T (GlobalParam_T &p);
	~GlobalParam_T (void);
	void setLamta (double *pLamta);
	void setDt (double dt);
	void setImageNum (int imageNum);
	void setMinObjVal (double eps);
	void setGapWidth (int gapWidth);
	void setWidth (int width);
	void setHeight (int height);
	void setMaxIter (int M) { maxIter = M; }
	void setSelfInit (bool s) { selfInit = s; }
	double *getLamta (void);
	double getDt (void);
	int getALength (void);
	int getImageNum (void);
	double getMinObjVal (void);
	int getGapWidth (void);
	int getWidth (void);
	int getHeight (void);
	int getMaxIter (void) { return maxIter; }
	bool getSelfInit (void) { return selfInit; } // modified from int to bool
private:
	double i_pLamta[INV_NUM];
	double i_dt;
	int i_aLength;
	int i_imageNum;
	double i_minObjVal;
	int i_gapWidth;
	int i_width;    // image width
	int i_height;   // image height
	int maxIter;
	bool selfInit;  // whether self initial {a}
};