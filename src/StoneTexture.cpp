
#include "StoneTexture.h"
#include "Perlin.h"


StoneTexture::StoneTexture() {
	StoneTexture(100);
}

StoneTexture::StoneTexture(int numCells) {
	stR = 160.0f/255.0f;
	stG = 82.0f/255.0f;
	stB = 45.0f/255.0f;

	//base color for grout
	gR = 250.0f/255.0f;
	gG = 235.0f/255.0f;
	gB = 215.0f/255.0f;


	m_numCells = numCells;
	xPts = new int[m_numCells];
	yPts = new int[m_numCells];
    //set up random points
	for (int i = numCells-1; i >= 0; i--) {
		xPts[i] = rand() % TEX_SIZE;
		yPts[i] = rand() % TEX_SIZE;
	}
	//get min and max for normalization
	m_minD = FLT_MAX;
	m_maxD = 0;
	for (int w = TEX_SIZE - 1; w >= 0; w--) {
		for (int h = TEX_SIZE - 1; h >= 0; h--) {
			float min = FLT_MAX;
			float min2 = FLT_MAX;
			for (int i = m_numCells - 1; i >= 0; i --) {
				float dist = pow((float)xPts[i]-w, 2) + pow((float)yPts[i] - h, 2);
				if (dist < min2) {
					if (dist < min) {
						min2 = min;
						min = dist;
					} else {
						min2 = dist;
					}
				}
			}
			min = sqrt(min2) - sqrt(min);
			if (min < m_minD)
				m_minD = min;
			if (min > m_maxD)
				m_maxD = min;
		}
	}
}

StoneTexture::~StoneTexture() {
	delete xPts;
	delete yPts;
}

Vector3 StoneTexture::getLookup3(float u, float v) {
	u = u - float(int(u));
	v = v - float(int(v));
	if (u < 0.0f)
		u = u + 1.0f;
	if (v < 0.0f)
		v = v + 1.0f;
	//textures start with v = 0 at the top, so reverse our v:
	v = 1.0f - v;

	float min = FLT_MAX;
	float min2 = FLT_MAX;
	int cellNum = 0;
	//find the F2 - F1 value
	for (int i = m_numCells - 1; i >= 0; i --) {
		float dist = pow((float)xPts[i]-u*TEX_SIZE, 2) + pow((float)yPts[i] - v*TEX_SIZE, 2);
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
	min = sqrt(min2) - sqrt(min);
	//stone mask is for the grout on the stone edges
	float stoneMask = (min - m_minD) / (m_maxD - m_minD);
	if (stoneMask > .05f) 
		stoneMask = 1.0f;
	else 
		stoneMask = 0.0f;

	//calculate the stone color
	float cn = .5*PerlinNoise::noise(255*float(cellNum)/m_numCells, 1.0f, 1.0f);
	float grout = .5 + .5*PerlinNoise::noise(255*u, 255*v, 1.0f);
	float stoneNoise = .05*PerlinNoise::noise(64*u, 64*v, 1.0f);
	Vector3 stoneColor;
	stoneColor.x = stoneNoise+(stR+cn) * stoneMask + gR*(1.0f - stoneMask)*grout;
	stoneColor.y = stoneNoise+(stG+cn*.2) * stoneMask + gG*(1.0f - stoneMask)*grout;
	stoneColor.z = stoneNoise+(stB+cn*.1) * stoneMask + gB*(1.0f - stoneMask)*grout;
	return stoneColor;
}

Vector4 StoneTexture::getLookup(float u, float v) {
	Vector3 temp = getLookup3(u,v);
	return Vector4(temp.x, temp.y, temp.z, 1.0f);
}