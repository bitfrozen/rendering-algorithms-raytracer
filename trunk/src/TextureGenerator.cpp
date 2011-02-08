#include <math.h>
#include "TextureGenerator.h"
#include "RawImage.h"
#include "Texture.h"


Texture* TextureGenerator::generateCellTexture(int width, int height, int numCells) {
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
		rawData[i] = pow(rawData[i], 2);
		//if (rawData[i] < .65f) rawData[i] = 0;
		//rawData[i] = pow(rawData[i], 2);
	}

	//create the new texture object
	RawImage* ri = new RawImage(width, height, rawData, GRAYSCALE);
	Texture* tex = new Texture(ri);

	
	//clean up
	delete xPts;
	delete yPts;
	return tex;
}


Texture* TextureGenerator::generateStoneTexture(int width, int height, int numCells) {
	float* stoneMask = new float[width*height];
	int (*xPts) = new int[numCells];
	int (*yPts) = new int[numCells];
	int* stoneColor = new int[width*height];
	float* rawData = new float[3*width*height];

	//base color for stone
	float stR = 160.0f/255.0f;
	float stG = 82.0f/255.0f;
	float stB = 45.0f/255.0f;

	//base color for grout
	float gR = 250.0f/255.0f;
	float gG = 235.0f/255.0f;
	float gB = 215.0f/255.0f;


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
			float min2 = FLT_MAX;
			int cellNum = 0;
			for (int i = numCells - 1; i >= 0; i --) {
				float dist = pow((float)xPts[i]-w, 2) + pow((float)yPts[i] - h, 2);
				if (dist < min2) {
					if (dist < min) {
						min2 = min;
						min = dist;
						cellNum = i;
					} else {
						min2 = dist;

					}

				}
			}
			//min = sqrt(min);
			min = sqrt(min2) - sqrt(min);
			//min = min2 - min;
			stoneMask[w*height + h] = min;
			stoneColor[w*height + h] = cellNum + 1;
			if (min < minDist)
				minDist = min;
			if (min > 200000){
				maxDist = 200000;
			}
			if (min > maxDist)
				maxDist = min;
		}
	}

	//normalize and calculate color
	for (int i = width*height - 1; i >= 0; i--) {
		int w = i/height;
		int h = i - w*height;
		//make the grout screen
		stoneMask[i] = (stoneMask[i] - minDist)/(maxDist - minDist);
		if (stoneMask[i] > .05f) 
			stoneMask[i] = 1.0f;
		else 
			stoneMask[i] = 0.0f;

	    //make the color of the stone:
		float cn = .25*PerlinNoise::noise(255*float(stoneColor[i])/numCells, 1.0f, 1.0f);
		float grout = .5 + .5*PerlinNoise::noise(255*float(w)/width, 255*float(h)/height, 1.0f);
		float stoneNoise = .05*PerlinNoise::noise(64*float(w)/width, 64*float(h)/height, 1.0f);
		rawData[i*3] = stoneNoise+(stR+cn) * stoneMask[i] + gR*(1.0f - stoneMask[i])*grout; 
		rawData[i*3+1] = stoneNoise+(stG+cn*.2) * stoneMask[i] + gG*(1.0f - stoneMask[i])*grout; 
		rawData[i*3+2] = stoneNoise+(stB+cn*.1) * stoneMask[i] + gB*(1.0f - stoneMask[i])*grout; 	
	}

	//create the new texture object
	RawImage* ri = new RawImage(width, height, rawData, RGB);
	Texture* tex = new Texture(ri);
	
	//clean up
	delete xPts;
	delete yPts;
	delete stoneMask;
	delete stoneColor;
	return tex;
}