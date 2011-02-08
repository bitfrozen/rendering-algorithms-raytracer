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
    Vector3 Ld = Vector3(0.0f, 0.0f, 0.0f);						// Diffuse	
	Vector3 Lr = Vector3(0.0f, 0.0f, 0.0f);						// Reflection
	Vector3 Lt = Vector3(0.0f, 0.0f, 0.0f);						// Transmission
	Vector3 Ls = Vector3(0.0f, 0.0f, 0.0f);						// Specular
    
	Vector3 rayD		= Vector3(ray.dx, ray.dy, ray.dz);		// Ray direction
    Vector3 viewDir		= -rayD;								// d is a unit vector
	float vDotN			= dot(viewDir, hit.N);					// Find if the interpolated normal points back while
	float vDotGeoN		= dot(viewDir, hit.geoN);				// geometric normal is fine.
	bool nEqGeoN		= (vDotN*vDotGeoN >= 0.0);				// True if both are equal
	Vector3 theNormal	= nEqGeoN ? hit.N : hit.geoN;			// Use geometric normal if interpolated one points away
	vDotN				= nEqGeoN ? vDotN : vDotGeoN;			// from object

	Vector3 P = Vector3(hit.px, hit.py, hit.pz);				// Hit position
	Vector3 rVec = (rayD + 2*(vDotN)*theNormal);				// get the reflection vector
    
    const Lights *lightlist = scene.lights();
	    
    // loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
        PointLight* pLight = *lightIter;
    
        Vector3 l = pLight->position() - P;
        
        // the inverse-squared falloff
        float falloff = l.length2();
		float distance = sqrt(falloff);
        
        // normalize the light direction
        l /= distance;

		float nDotL = std::max(0.0f, dot(theNormal, l));
		Vector3 E = (nDotL * pLight->color() * pLight->wattage()) / (4*falloff*PI);		// Light irradiance
		Ray sRay = Ray(P, l, 1.0f, 0, IS_SHADOW_RAY);									// Create shadow ray
		HitInfo sHit;
		if (nDotL)
		{
			if (scene.trace(sHit, sRay, 0.001, distance))								// Check for shadows
			{
				E *= (sHit.t < distance) ? 0 : 1;
			}
		}	

		Ls += m_ks*m_specAmt*pow(dot(rVec, l), m_specExp);								// Calculate specular component
		Ld += E*(m_kd + Ls);															// Total = Light*(Diffuse+Spec)
    }

	bool flip;
	Ray::IORList IORHistory = ray.r_IOR;
	float outIOR;									// outIOR is the IOR of the medium the ray will go into.
	float inIOR = IORHistory.back();				// inIOR is the IOR of the medium the ray is currently in.
	if (vDotN < 0.0)								// If this is true, then we've hit a back-facing polygon:
	{												// we're going out of the current material.
		if (IORHistory.size() > 1)					// Pop the back of the IORHistory so we get the IOR right, be careful and make sure this is OK.
		{
			IORHistory.pop_back();
		}		
		flip = true;
		vDotN = -vDotN;								// Flip the vDotN so that fresnel and angle calculations are right
		outIOR = IORHistory.back();
	}
	else											// Normal. We're going (if possible) into a new material,
	{												// record its IOR in the ray's IORHistory
		flip = false;
		outIOR = m_ior;
	}
	float Rs = (m_reflectAmt > 0.0f) ? Material::fresnel(inIOR, outIOR, vDotN) : 0.0;	// Compute Fresnel's equations  http://www.bramz.net/data/writings/reflection_transmission.pdf
	float Ts = 1.0f-Rs;																	// Energy conservation
	int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);						// Find out how many bounces we've taken

	if (m_reflectAmt*Ts > 0.0f && bounces < 2)											// Do two bounces of reflection rays.
	{
		Ray rRay = Ray(P, rVec, ray.r_IOR, bounces+1, IS_REFLECT_RAY);
		HitInfo rHit;
		if (scene.trace(rHit, rRay, 0.001, MIRO_TMAX))
		{
			Vector3 reflection = rHit.material->shade(rRay, rHit, scene);				
			Lr += m_ks*m_reflectAmt*Rs*reflection;
		}
		else Lr += m_ks*m_reflectAmt*Rs*g_scene->getBGColor();
	}
	else Lr += m_ks*m_reflectAmt*Rs*g_scene->getBGColor();

	if (m_refractAmt*Ts > 0.0f && bounces < 4)											// Do two bounces of refraction rays (2 bounces have already been used by reflection).
	{
		float snellsQ = inIOR / outIOR;
		float sqrtPart = sqrt(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN));
		Vector3 tVec;
		if (flip)
		{
			tVec = snellsQ*rayD - theNormal*(snellsQ*vDotN - sqrtPart);					// Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf
		} 
		else
		{
			tVec = snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart);
		}		
		IORHistory.push_back(outIOR);
		Ray tRay = Ray(P, tVec, IORHistory, bounces+1, IS_REFRACT_RAY);					// Make a new refraction ray.
		HitInfo tHit;
		if (scene.trace(tHit, tRay, 0.001, MIRO_TMAX))
		{
			Vector3 refraction = tHit.material->shade(tRay, tHit, scene);				
			Lt += m_ks*m_refractAmt*Ts*refraction;
		}
		else Lt += m_ks*m_refractAmt*Ts*g_scene->getBGColor();
	}
	else Lt += m_ks*m_refractAmt*Ts*g_scene->getBGColor();
    
    Ld += m_ka;
    
    return Ld*(1.0-Rs*m_reflectAmt-Ts*m_refractAmt)+Lr+Lt;
}
