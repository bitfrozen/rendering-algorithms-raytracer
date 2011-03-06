#ifndef CSE168_DOMELIGHT_H_INCLUDED
#define CSE168_DOMELIGHT_H_INCLUDED

#include "Light.h"
#include "Texture.h"
#include "Scene.h"
#include <algorithm>

// Based on implementation details from PBRT
struct Distribution1D {
	Distribution1D(float *f, int n)
	{
		func = new float[n];
		cdf = new float[n+1];
		count = n;
		memcpy(func, f, n*sizeof(float));
		computeStep1dCDF(func, n, &funcInt, cdf);
		invFuncInt = 1.f / funcInt;
		invCount = 1.f / count;
	}
	void computeStep1dCDF(float *f, int nSteps, float *c, float *cdf)
	{
		int i;
		cdf[0] = 0.f;
		for (i = 1; i < nSteps+1; ++i)
			cdf[i] = cdf[i-1] + f[i-1] / nSteps;
		*c = cdf[nSteps];
		for (i = 1; i < nSteps+1; ++i)
			cdf[i] /= *c;		
	}
	float sample(float u, float *pdf)
	{
		float *ptr = std::lower_bound(cdf, cdf+count+1, u);
		int offset = (int) (ptr-cdf-1);
		u = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
		*pdf = func[offset] * invFuncInt;
		return offset + u;
	}
	float *func, *cdf;
	float funcInt, invFuncInt, invCount;
	int count;
};

// Implements a direct sampled arbitrary parallelogram light
class DomeLight : public Light
{
public:
	DomeLight() : Light() {m_lightMap = NULL; m_Color = Vector3(0); m_Gain = 1.f; setType(DOME_LIGHT); cosTableU = sinTableU = cosTableV = sinTableV = NULL; uDistrib = NULL; vDistribs = NULL;}

	void setTexture(Texture* t);
	void setColor(const Vector3& v) {m_Color = v;}
	void setPower(float f) {m_Gain = f;}

	const Vector3 sampleLight(const unsigned int threadID, const Vector3 &from, const Vector3 &normal, const float time, const Scene &scene, const Vector3 &rVec, float &outSpec, bool isSecondary = false) const;

protected:
	Texture* m_lightMap;
	float m_Gain;
	Vector3 m_Color;
	Distribution1D *uDistrib, **vDistribs;
	float *cosTableU, *cosTableV, *sinTableU, *sinTableV;
};

#endif