#ifndef CSE_168_HDRLOADER
#define CSE_168_HDRLOADER
/***********************************************************************************
	Created:	17:9:2002
	FileName: 	hdrloader.h
	Author:		Igor Kravtchenko
	
	Info:		Load HDR image and convert to a set of float32 RGB triplet.
************************************************************************************/
/*
class HDRLoaderResult {
public:
	int width, height;
	// each pixel takes 3 float32, each component can be of any value...
	float *cols;
};*/

#include "RawImage.h"

class HDRLoader {
public:
	static bool load(const char *fileName, float* &fdata, int &width, int &height);
};


#endif