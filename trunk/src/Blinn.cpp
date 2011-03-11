#include "Blinn.h"
#include "Ray.h"
#include "Scene.h"
#include <time.h>
#include <math.h>
#include "Miro.h"
#include "SSE.h"
#include "TriangleMesh.h"
#include "RectangleLight.h"
#include "PointLight.h"
#include "Camera.h"

using namespace std;

Blinn::Blinn(const Vector3 & kd, const Vector3 & ka,
			 const Vector3 & ks, const Vector3 & kt,
			 float ior, float specExp, float specAmt,
			 float reflectAmt, float refractAmt, float specGloss,
			 float emittedPower, const Vector3& emittedColor) :
    m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt),
	m_specExp(specExp), m_specAmt(specAmt),
	m_reflectAmt(reflectAmt), m_specGloss(specGloss),
	m_lightEmitted(emittedPower), m_Le(emittedColor)
{	
	m_ior[0] = ior;
	m_ior[1] = ior;
	m_ior[2] = ior;
	m_lightEmitted = 0.0f;
	m_Le = Vector3(0.0f);
	m_refractAmt = refractAmt;
	m_translucency = 0.0f;
	m_disperse = false;
}

Blinn::~Blinn()
{
}

const Vector3 Blinn::calculatePathTracing(const unsigned int threadID, const Ray& ray, const HitInfo &hit, const Vector3& P, const Vector3& theNormal, const Scene& scene, const Vector3& diffuseColor) const
{
	Vector3 out = 0;
	Vector3 randD, envColor;
	Ray randRay(threadID);
	HitInfo newHit;

	// If we get to an emitter (mesh light), we return its intensity
	if (m_lightEmitted > 0.0f || (m_Le.x + m_Le.y + m_Le.z) > 0.0f)
	{
		out += m_lightEmitted * m_Le;
		return out;
	}

	// Compute GI component using path tracing.
	int bounces   = GET_BOUNCES(ray.bounces_flags);
	int giBounces = GET_GI_BOUNCES(ray.bounces_flags);

	if (giBounces < scene.m_maxBounces-1) 
	{
		getCosineDistributedSamples(threadID, theNormal, randD);

		randRay.set(threadID, P, randD, ray.time, ray.r_IOR(), bounces, giBounces+1, IS_PRIMARY_RAY);
		newHit.t = MIRO_TMAX;

		// Trace the random ray; Do the next bounce
		if (scene.trace(threadID, newHit, randRay, epsilon))
		{
			out += diffuseColor * newHit.obj->m_material->shade(threadID, randRay, newHit, scene, true);
		}
		// If we don't hit anything, add contribution from environment
		else if (m_sampleEnv && g_scene->sampleEnv())
		{
			out += diffuseColor * getEnvironmentColor(randD, scene);
		}
	}
	// We directly sample the light at the last bounce...
	else
	{
		// Directly sample lights
		const Lights *lightList = scene.lights();
		for (int i = 0; i < lightList->size(); i++)
		{
			float   lightSpec  = 0;
			Vector3 lightPower = (*lightList)[i]->sampleLight(threadID, P, theNormal, ray.time, scene, Vector3(0), lightSpec, true);
		
			out += lightPower * diffuseColor;										// Calculate Diffuse component
		}
	}
	return out;
}

const Vector3 Blinn::shade(const unsigned int threadID, const Ray& ray, const HitInfo &hit, const Scene& scene, bool isSecondary) const
{
	Vector3 Ld		= 0.0f;										// Diffuse	
	Vector3 Lr		= 0.0f;										// Reflection
	Vector3 Lt		= 0.0f;										// Transmission
	Vector3 Ls		= 0.0f;										// Specular
	Vector3 rayD	= Vector3(ray.d[0],ray.d[1],ray.d[2]);		// Ray direction
	Vector3 viewDir	= -rayD;									// View direction
	float u, v;
	Vector3 N, geoN, T, BT;
	Vector3 diffuseColor = m_kd;
	float localSpecExp = m_specExp;
	float localSpecAmt = m_specAmt;
	float localReflectAmt = m_reflectAmt;
	float localRefractAmt = m_refractAmt;
	float alpha = 1.0f;

	hit.getAllInfos(N, geoN, T, BT, u, v);
	Vector3 P = ray.getPoint(hit.t);
	int bounces   = GET_BOUNCES(ray.bounces_flags);
	int giBounces = GET_GI_BOUNCES(ray.bounces_flags);
	Vector3 translucency(0.0f);

	if (m_colorMap) 
	{
		Vector4 texCol = m_colorMap->getLookup(u, v);
		diffuseColor = Vector3(texCol.x, texCol.y, texCol.z);
	}

	if (m_normalMap)
	{
		Vector4 texN = m_normalMap->getLookup(u, v);
		N = (texN.x*T + texN.y*BT + texN.z*N);
	}

	if (m_specularMap)
	{
		Vector4 texS = m_specularMap->getLookup(u, v);
		localSpecAmt = (texS.x+texS.y+texS.z) * 0.3333333f * localSpecAmt;
	}

	if (m_reflectMap)
	{
		Vector4 texR = m_reflectMap->getLookup(u, v);
		localReflectAmt = (texR.x+texR.y+texR.z) * 0.3333333f * localReflectAmt;
	}

	if (m_refractMap)
	{
		Vector4 texR = m_refractMap->getLookup(u, v);
		localRefractAmt = (texR.x+texR.y+texR.z) * 0.3333333f * localRefractAmt;
	}

	bool flip = false;
	float vDotN			= dot(viewDir, N);						// Find if the interpolated normal points back while
	float vDotGeoN		= dot(viewDir, geoN);					// geometric normal is fine.
	bool nEqGeoN		= (vDotN*vDotGeoN >= 0.0);				// True if both are equal
	Vector3 theNormal	= nEqGeoN ? N : geoN;					// Use geometric normal if interpolated one points away
	vDotN				= nEqGeoN ? vDotN : vDotGeoN;			// from object
	if (vDotN < 0.0)
	{
		flip = true;
		vDotN = -vDotN;
		theNormal = -theNormal;
	}
	
	// get the reflection vector (for specular hi-light)
	Vector3 rVec = (rayD + 2*(vDotN)*theNormal);

	Vector3 randD;																// If we have glossy reflections, randomize the reflection vector.
	if (m_specGloss < 1.0)
	{
		getCosineDistributedSamples(threadID, theNormal, randD);
		rVec = (m_specGloss*rVec + (1-m_specGloss)*randD).normalized();			// get the (randomized) reflection vector
	}

	float outIOR[3];									// outIOR is the IOR of the medium the ray will go into.
	float inIOR = ray.r_IOR();						// inIOR is the IOR of the medium the ray is currently in.
	if (m_disperse && !(ray.bounces_flags & IS_REFRACT_RAY)){
		outIOR[0] = m_ior[0];
		outIOR[1] = m_ior[1];
		outIOR[2] = m_ior[2];
	}
	else
	{
		if (flip)										// If this is true, then we've hit a back-facing polygon:
		{												// we're going out of the current material.
			ray.r_IOR.pop();							// Pop the back of the IORHistory so we get the IOR right.
			outIOR[0] = ray.r_IOR();
		}
		else											// Normal. We're going (if possible) into a new material,
		{												// record its IOR in the ray's IORHistory
			outIOR[0] = m_ior[1];
		}
	}

	// Calculate reflectance and transmission coefficients
	float Rs = 0; float Ts = 0;
	if (m_reflectAmt > 0.0 || m_refractAmt > 0.0)
	{
		Rs = fresnel(inIOR, outIOR[0], vDotN);				// Compute Fresnel's equations http://www.bramz.net/data/writings/reflection_transmission.pdf
		Ts = 1.0f-Rs;										// Energy conservation
	}

	float rrFloat = Scene::getRand(threadID);
	float rrWeight = (1.0f - Rs*localReflectAmt - Ts*localRefractAmt);
	float rrWeightRecip = (rrWeight > 0.f) ? 1.f / rrWeight : 1.f;
	float rrWeightRecipSpec = (1.f-rrWeight > 0.f) ? 1.f / (1.f-rrWeight) : 1.f;

	// Russian roulette sampling of direct vs specular illumination
	// This is weighted by the respective diffuse/specular balance of
	// the material at the shaded point
	
	if (rrFloat <= rrWeight)
	{
		// Get the Indirect Illumination
		if (scene.m_pathTrace)
		{	
			Ld += calculatePathTracing(threadID, ray, hit, P, theNormal, scene, diffuseColor);// * 2.f;
		}

		// Directly sample lights
		const Lights *lightList = scene.lights();
		for (int i = 0; i < lightList->size(); i++)
		{
			float   lightSpec  = 0;
			Vector3 lightPower = (*lightList)[i]->sampleLight(threadID, P, theNormal, ray.time, scene, rVec, lightSpec, isSecondary);

			Ls += lightPower * m_ks * localSpecAmt * pow(lightSpec, localSpecExp);// * 2.f;		// Calculate specular component
			Ld += lightPower * diffuseColor;// * 2.f;										// Calculate Diffuse component
		}

		if (m_translucency > 0.01f) 
		{
			//sample light to get translucent color, or should we shot more rays?
			Vector3 lightTotal(0.0f);
			const Lights *lightList = scene.lights();
			for (int i = 0; i < lightList->size(); i++)
			{
				float   lightSpec  = 0;
				Vector3 lightPower = (*lightList)[i]->sampleLight(threadID, P, -theNormal, .001f, scene, rVec, lightSpec, isSecondary);
				lightTotal += lightPower;
			}
			//update the light with the new color
			translucency += m_translucency * lightTotal * diffuseColor;
		}
	}	
	else
	{
		bool doEnv  = true;

		rrFloat = Scene::getRand(threadID);

		// Russian Roulette sampling for reflection / refraction.
		// This helps prevent ray number explosion
 		if (rrFloat < localReflectAmt*Rs)
 		{
			if (localReflectAmt*Rs > 0.0f && bounces < 5)
			{
				Ray rRay = Ray(threadID, P, rVec, ray.time, ray.r_IOR, bounces+1, giBounces, IS_REFLECT_RAY);
				HitInfo newHit;
				newHit.t = MIRO_TMAX;
				if (scene.trace(threadID, newHit, rRay, epsilon))
				{
					Vector3 reflection = newHit.obj->m_material->shade(threadID, rRay, newHit, scene);				
					Lr                += m_ks*reflection;
					doEnv              = false;
				}
			}
			if (localReflectAmt*Rs > 0.0f && doEnv)
			{
				Lr += m_ks * getEnvironmentColor(rVec, scene);
				if (Lr.average() > 10000.f)
				{
					int tmp = 0;
				}
			}
 		} 
 		else
 		{
			if (localRefractAmt*Ts > 0.0f)
			{  
				doEnv          = true;
				Vector3 tVec;
				if (m_disperse && !(ray.bounces_flags & IS_REFRACT_RAY)) 
				{
					for (int i = 0; i < 3; i++) {
						//0 = r, 1 = g, 2 = b
						float snellsQ  = inIOR / outIOR[i];
						float sqrtPart = max(0.0f, sqrtf(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN)));
						tVec   = (snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart)).normalized();	// Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf

						if (bounces < 5)																		// Do two bounces of refraction rays (2 bounces have already been used by reflection).
						{		
							ray.r_IOR.push(outIOR[i]);//shot out 3+ rays, dont split up if IS_REFLECT_RAY is set
							Ray tRay = Ray(threadID, P, tVec, ray.time, ray.r_IOR, bounces+1, giBounces, IS_REFRACT_RAY);			// Make a new refraction ray.
							HitInfo newHit;
							newHit.t = MIRO_TMAX;
							if (scene.trace(threadID, newHit, tRay, epsilon))
							{
								Vector3 refraction = newHit.obj->m_material->shade(threadID, tRay, newHit, scene);
								Vector3 mask(0.0f);
								mask[i] = 1.0f;
								refraction         = refraction*mask;
								Lt                += m_ks*refraction;
								doEnv              = false;
							}
							ray.r_IOR.pop();
						}
					}
				} 
				else 
				{
					//no dispersion
					float snellsQ  = inIOR / outIOR[0];
					float sqrtPart = max(0.0f, sqrtf(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN)));
					tVec   = (snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart)).normalized();	// Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf

					if (bounces < 5)																		// Do two bounces of refraction rays (2 bounces have already been used by reflection).
					{		
						ray.r_IOR.push(outIOR[0]);
						Ray tRay = Ray(threadID, P, tVec, ray.time, ray.r_IOR, bounces+1, giBounces, IS_REFRACT_RAY);			// Make a new refraction ray.
						HitInfo newHit;
						newHit.t = MIRO_TMAX;
						if (scene.trace(threadID, newHit, tRay, epsilon))
						{
							Vector3 refraction = newHit.obj->m_material->shade(threadID, tRay, newHit, scene);				
							Lt                += m_ks*refraction;
							doEnv              = false;
						}
						ray.r_IOR.pop();
					}
				}

				if (doEnv)
				{
					Lt += m_ks * getEnvironmentColor(tVec, scene);
				}
			}
		}
	}

	Ld += m_ka;

	return (Ld + Ls + translucency)*rrWeightRecip + (Lr + Lt)*rrWeightRecipSpec + m_Le; // Make sure we keep energy constant
}