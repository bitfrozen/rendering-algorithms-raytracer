#include "Scene.h"
#include "Camera.h"
#include "PointLight.h"
#include "SSE.h"

using namespace std;

const Vector3 PointLight::sampleLight(const unsigned int threadID, const Vector3 &from, const Vector3 &normal, const float time, const Scene &scene, const Vector3 &rVec, float &outSpec) const
{
	Ray sampleRay(threadID);
	HitInfo sampleHit;
	Vector3 L		= 0.0f;
	Vector3 E		= 0.0f;
	float attenuate = 1.0f;
	ALIGN_SSE float falloff = 0.0;

	// Get a vector into the light
	L = m_position - from;
	float nDotL = dot(normal, L);

	if (nDotL > 0.0f)															// Only do work if light can be hit
	{    
		// the inverse-squared falloff
		falloff = L.length2();
		ALIGN_SSE float distanceRecip, distance;
#ifndef NO_SSE
		fastrsqrtss(setSSE(falloff), distanceRecip);
		recipss(setSSE(falloff), falloff);
		recipss(setSSE(distanceRecip), distance);
#else
		distance      = sqrtf(falloff);
		distanceRecip = 1.0f / distance;
		falloff       = 1.0f / falloff;
#endif
		L     *= distanceRecip;
		nDotL *= distanceRecip;

		if (m_castShadows)
		{
			sampleHit.t = distance;
			if (m_fastShadows)
			{
				sampleRay.set(threadID, from, L, time, 1.001f, 0, IS_SHADOW_RAY);			// Create shadow ray
				if (scene.trace(threadID, sampleHit, sampleRay, 0.001))				// Quick method, returns any hit
				{
					attenuate = 0.0f;
				}
			}
			else																// Full method, accounts for transparency effects
			{
				sampleRay.set(threadID, from, L, time, 1.001f, 0, IS_PRIMARY_RAY);	// Create primary ray so we trace properly
				while (sampleHit.t < distance)
				{
					if (scene.trace(threadID, sampleHit, sampleRay, 0.001))				
					{
						Vector3 hitN; sampleHit.getInterpolatedNormal(hitN);
						float nDL = dot(hitN, -L);
						if (nDL > 0.0)											// Only attenuate on incoming direction
						{
							attenuate *= sampleHit.obj->m_material->refractAmt();
						}
						Vector3 newPoint = Vector3(sampleRay.o[0], sampleRay.o[1], sampleRay.o[2]) + sampleHit.t * L;
						sampleRay.set(threadID, newPoint, L, time, 1.001f, IS_PRIMARY_RAY);
					}
					else
					{
						sampleHit.t = MIRO_TMAX;
					}
				}
			}
		}
		attenuate *= nDotL;														// Take the cosine term into account
	}
	else
	{
		outSpec = 0;
		return Vector3(0);
	}

	outSpec = max(0.f, dot(rVec, L)) * attenuate;
	return E = m_power * falloff *_1_4PI * attenuate;							// Light irradiance for this sample;
}