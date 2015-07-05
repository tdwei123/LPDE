
// IMPORTANT NOTES: function setArea must be called before using an instance of DifImage_T class

#pragma once

class DifImage_T {
public:
	DifImage_T(int width = 0, int height = 0);
	~DifImage_T(void);
	void setArea(int width, int height);
	void setVal(double *pImageVal);
	void setVal(DifImage_T &img);
	double getVal(int i, int j, int x, int y);
	double getVal(int i, int j, int pos);
	double *getVal(int i, int j);
	int getWidth() { return i_width;}
	int getHeight() { return i_height;}
private:
	double *i_pVal;
	int i_width, i_height;
};