#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include "Miro.h"
#include "Vector3.h"
#include "Texture.h"

using namespace std;

class Material
{
public:
    Material();
    virtual ~Material();

    virtual void preCalc() {}
    
    const virtual Vector3 shade(const unsigned int threadID, const Ray& ray, const HitInfo& hit, const Scene& scene) const = 0;
	void setEnvMap(Texture* map)					{m_envMap = map;}
	void setEnvExposure(float exp)					{m_envExposure = exp;}
	const float refractAmt() const					{return m_refractAmt;}

	void setRefractAmt(const float refractAmt)		{m_refractAmt = refractAmt;}
	const Vector3 getEnvironmentColor(const Vector3& direction, const Scene& scene) const;

	Texture* m_texture;
	Texture* m_envMap;
	float m_envExposure;

	__forceinline static float fresnel(const float n1, const float n2, const float cosThetaI) 
	{
#ifndef USE_SCHLICK
		const float n1CosTh = n1*cosThetaI;
		const float n1_n2SinTh = n1*sin(acosf(cosThetaI))/n2;
		const float n2CosTh = n2*max(0.0f, sqrt(1.0f - n1_n2SinTh*n1_n2SinTh));
		const float Rs = (n1CosTh - n2CosTh) / (n1CosTh + n2CosTh);
		return Rs*Rs;
#else
		float r0 = (n1 - n2) / (n1 + n2);
		r0 *= r0;
		if (n1 > n2)
		{
			const float n = n1 / n2;
			const float sinT2 = n*n*(1.0 - cosThetaI*cosThetaI);
			if (sinT2 > 1.0) return 1.0; // TIR
			cosThetaI = sqrt(1.0 - sinT2);			
		}
		const float x = 1.0 - cosThetaI;
		const float x2 = x*x;
		return r0 + (1.0 - r0) * x2*x2*x;
#endif
	}

	static void getCosineDistributedSamples(const Vector3 &N, Vector3 &out);
protected:
	float m_refractAmt;				// Refraction amount. This weights the refraction amount prescribed by the fresnel approximation.
};

#endif // CSE168_MATERIAL_H_INCLUDED
