#include <math.h>
#include "TextureGenerator.h"
#include "RawImage.h"
#include "Texture.h"


Texture* const TextureGenerator::generateCellTexture(int width, int height, int numCells, Texture* noise) {
	float* rawData = new float[width*height];
	int (*xPts) = new int[numCells];
	int (*yPts) = new int[numCells];
	//int xPts[numCells];
	//int yPts[numCells];
	float minDist = FLT_MAX;
	float maxDist = 0;
	//initialize the random points on the grid
	for (int i = numCells-1; i >= 0; i--) {
	  xPts[i] = rand() % width;
	  yPts[i] = rand() % height;
	}
	//get distances to points
	for (int w = width - 1; w >= 0; w--) {
		for (int h = height - 1; h >= 0; h--) {
			float min = FLT_MAX;
			for (int i = numCells - 1; i >= 0; i --) {
				float dist = pow((float)xPts[i]-w, 2) + pow((float)yPts[i] - h, 2);
				if (dist < min) 
					min = dist;
			}
			min = sqrt(min);
			rawData[h*width + w] = min;
			if (min < minDist)
				minDist = min;
			if (min > maxDist)
				maxDist = min;
		}
	}

	//normalize
	for (int i = width*height - 1; i >= 0; i--) {
		rawData[i] = (rawData[i] - minDist)/(maxDist - minDist);
	}

	//create the new texture object
	RawImage* ri = new RawImage(width, height, rawData, GRAYSCALE);
	Texture* tex = new Texture(ri);
	
	//clean up
	delete xPts;
	delete yPts;
	return tex;
}
