#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include "Miro.h"
#include "Vector3.h"
#include "Texture.h"

class Material
{
public:
    Material();
    virtual ~Material();

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;
	void setEnvMap(Texture* map)	{m_envMap = map;}
	void setEnvExposure(float exp)	{m_envExposure = exp;}
	Texture* m_texture;
	Texture* m_envMap;
	float m_envExposure;

	static float fresnel(float n1, float n2, float cosThetaI) {		// Fn to compute fresnel coefficients.
		if (!use_Schlick)
		{
			float n1CosTh = n1*cosThetaI;
			float n1_n2SinTh = n1*sin(acosf(cosThetaI))/n2;
			float n2CosTh = n2*std::max(0.0f, sqrt(1.0f - n1_n2SinTh*n1_n2SinTh));
			float Rs = (n1CosTh - n2CosTh) / (n1CosTh + n2CosTh);
			return Rs*Rs;
		}
		float Rs = (n1-n2)/(n1+n2);
		Rs = Rs*Rs;
		return Rs + (1-Rs)*pow(1-cosThetaI, 5);
	}
};

#endif // CSE168_MATERIAL_H_INCLUDED
