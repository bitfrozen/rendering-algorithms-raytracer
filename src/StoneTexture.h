#ifndef CSE_168_STONETEXTURE
#define CSE_168_STONETEXTURE

#include "Vector3.h"
#include "Texture.h"

class StoneTexture : public Texture {
public:
	StoneTexture();
	StoneTexture(int numCells);
	~StoneTexture();
	virtual Vector3 getLookup3(float u, float v);
	virtual Vector4 getLookup(float u, float v);

private:
	int m_numCells;
	int* xPts;
	int* yPts;

	float m_minD, m_maxD;

	//base color for stone
	float stR, stG, stB;

	//base color for grout
	float gR, gG, gB;
	
	static const int TEX_SIZE = 256;
};


#endif