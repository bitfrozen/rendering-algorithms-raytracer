#include "Blinn.h"
#include "Ray.h"
#include "Scene.h"
#include <math.h>

Blinn::Blinn(const Vector3 & kd, const Vector3 & ka,
			 const Vector3 & ks, const Vector3 & kt,
			 float ior, float specExp, float specAmt,
			 float reflectAmt, float refractAmt) :
    m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt), m_ior(ior), 
	m_specExp(specExp), m_specAmt(specAmt),
	m_reflectAmt(reflectAmt), m_refractAmt(refractAmt)
{
}

Blinn::~Blinn()
{
}

Vector3
Blinn::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
    Vector3 Ld = Vector3(0.0f, 0.0f, 0.0f);
	Vector3 Lr = Vector3(0.0f, 0.0f, 0.0f);
	Vector3 Lt = Vector3(0.0f, 0.0f, 0.0f);
	Vector3 Ls = Vector3(0.0f, 0.0f, 0.0f);
    
    Vector3 viewDir = -ray.d; // d is a unit vector
	Vector3 rVec = (ray.d - 2*(dot(ray.d, hit.N))*hit.N);  // get the reflection vector
    
    const Lights *lightlist = scene.lights();
	    
    // loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;
    
        Vector3 l = pLight->position() - hit.P;
        
        // the inverse-squared falloff
        float falloff = l.length2();
        
        // normalize the light direction
        l /= sqrt(falloff);

		float nDotL = std::max(0.0f, dot(hit.N, l));
		Vector3 E = (nDotL * pLight->color() * pLight->wattage()) / (4*falloff*PI);
		Ray sRay = Ray(hit.P+0.002*l, l);
		HitInfo sHit;
		if (nDotL)
		{
			if (scene.trace(sHit, sRay))
			{
				E *= (sHit.t < sqrt(falloff)) ? 0 : 1;
			}
		}	

		Ls += m_ks*m_specAmt*pow(dot(rVec, l), m_specExp) * E;        
		Ld += m_kd*E;
    }

	float nDotV = dot(hit.N, viewDir); 
	float fresnel = std::min(m_reflectAmt + (1.0f - m_reflectAmt)*pow(1.0f - nDotV, 4.0f), 1.0f);

	/*if (m_reflectAmt > 0.0 && ray.bounces < 2)
	{
		Ray rRay = Ray(hit.P+0.002*rVec, rVec, ray.bounces+1);
		HitInfo rHit;
		if (scene.trace(rHit, rRay) && rHit.t < MIRO_TMAX)
		{
			Vector3 reflection = rHit.material->shade(rRay, rHit, scene);				
			Lr += m_ks*fresnel*reflection;
		}
		else Lr += m_ks*fresnel*scene.envColor();
	}*/
	/*
	if (m_refractAmt > 0.0 && ray.bounces < 4)
	{
		float inv_IOR = 1.0f / m_ior;
		float cosTh2 = sqrt(1.0f - (inv_IOR*inv_IOR) * (1.0f-nDotV*nDotV));
		Vector3 tVec = inv_IOR*ray.d + hit.N * ((nDotV > 0.0) ? (inv_IOR*nDotV-cosTh2) : (inv_IOR*nDotV+cosTh2));
		Ray tRay = Ray(hit.P+0.002*tVec, tVec, ray.bounces+1);
		HitInfo tHit;
		if (scene.trace(tHit, tRay) && tHit.t < MIRO_TMAX)
		{
			Vector3 reflection = tHit.material->shade(tRay, tHit, scene);				
			Ld = reflection;//m_ks*(1.f-fresnel)*reflection;
		}
		//else Ld += m_ks*(1.f-fresnel)*scene.envColor();
	}*/
    
    // add the ambient component
    //Ld += m_ka;
    
    return Ld;//+Lr+Ls;
}
