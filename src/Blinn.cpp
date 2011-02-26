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

using namespace std;

Blinn::Blinn(const Vector3 & kd, const Vector3 & ka,
			 const Vector3 & ks, const Vector3 & kt,
			 float ior, float specExp, float specAmt,
			 float reflectAmt, float refractAmt, float specGloss,
			 float emittedPower, const Vector3& emittedColor) :
    m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt), m_ior(ior), 
	m_specExp(specExp), m_specAmt(specAmt),
	m_reflectAmt(reflectAmt), m_specGloss(specGloss),
	m_lightEmitted(emittedPower), m_Le(emittedColor)
{	
	m_lightEmitted = 0.0f;
	m_Le = Vector3(0.0f);
	m_refractAmt = refractAmt;
}

Blinn::~Blinn()
{
}

const Vector3 Blinn::calculatePathTracing(const Ray& ray, const HitInfo& hit, const Vector3& P, const Vector3& theNormal, const Scene& scene, const Vector3& diffuseColor) const
{
	Vector3 out = 0;
	Vector3 randD, envColor;
	Ray randRay;
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

	if (giBounces < scene.m_maxBounces) 
	{
		getCosineDistributedSamples(theNormal, randD);

		randRay.set(P, randD, ray.r_IOR(), bounces, giBounces+1, IS_PRIMARY_RAY);
		newHit.t = MIRO_TMAX;

		// Trace the random ray; Do the next bounce
		if (scene.trace(newHit, randRay, epsilon))
		{
			out += diffuseColor * newHit.obj->m_material->shade(randRay, newHit, scene);
		}
		// If we don't hit anything, add contribution from environment
		else
		{
			out += diffuseColor * getEnvironmentColor(randD, scene);
		}
	}	
	return out;
}

Vector3 Blinn::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
	Vector3 Ld		= 0.0f;										// Diffuse	
	Vector3 Lr		= 0.0f;										// Reflection
	Vector3 Lt		= 0.0f;										// Transmission
	Vector3 Ls		= 0.0f;										// Specular
	Vector3 rayD	= Vector3(ray.d[0],ray.d[1],ray.d[2]);		// Ray direction
	Vector3 viewDir	= -rayD;									// View direction
	float u, v;
	Vector3 N, geoN;
	Vector3 diffuseColor = m_kd;

	hit.getAllInfos(N, geoN, u, v);

	if (m_texture != NULL) 
	{
		Vector4 texCol = m_texture->getLookup(u, v);
		diffuseColor = Vector3(texCol.x, texCol.y, texCol.z);
	}

	Vector3 P = ray.getPoint(hit.t);

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
	
	// Get the Indirect Illumination
	if (scene.m_pathTrace)
	{	
		Ld += calculatePathTracing(ray, hit, P, theNormal, scene, diffuseColor);
	}
	
	// get the reflection vector (for specular hi-light)
	Vector3 rVec = (rayD + 2*(vDotN)*theNormal);

	Vector3 randD;																// If we have glossy reflections, randomize the reflection vector.
	getCosineDistributedSamples(theNormal, randD);
	if (m_specGloss < 1.0)
	{
		rVec = (m_specGloss*rVec + (1-m_specGloss)*randD).normalized();			// get the (randomized) reflection vector
	}

	// Directly sample lights
	const Lights *lightList = scene.lights();
	Lights::const_iterator lightIter;
	for (lightIter = lightList->begin(); lightIter != lightList->end(); lightIter++)
	{
		float   lightSpec  = 0;
		Vector3 lightPower = (*lightIter)->sampleLight(P, theNormal, scene, rVec, lightSpec);

		Ls += lightPower * m_ks * m_specAmt * pow(lightSpec, m_specExp);		// Calculate specular component				
		Ld += lightPower * diffuseColor;										// Calculate Diffuse component
	}		

	float outIOR;									// outIOR is the IOR of the medium the ray will go into.
	float inIOR = ray.r_IOR();						// inIOR is the IOR of the medium the ray is currently in.
	if (flip)										// If this is true, then we've hit a back-facing polygon:
	{												// we're going out of the current material.
		ray.r_IOR.pop();							// Pop the back of the IORHistory so we get the IOR right.
		outIOR = ray.r_IOR();
	}
	else											// Normal. We're going (if possible) into a new material,
	{												// record its IOR in the ray's IORHistory
		outIOR = m_ior;
	}

	// Calculate reflectance and transmission coefficients
	float Rs = 0; float Ts = 0;
	if (m_reflectAmt > 0.0 || m_refractAmt > 0.0)
	{
		Rs = fresnel(inIOR, outIOR, vDotN);					// Compute Fresnel's equations http://www.bramz.net/data/writings/reflection_transmission.pdf
		Ts = 1.0f-Rs;										// Energy conservation
	}
	
	int bounces   = GET_BOUNCES(ray.bounces_flags);
	int giBounces = GET_GI_BOUNCES(ray.bounces_flags);

	bool doEnv  = true;

	if (m_reflectAmt*Rs > 0.0f && bounces < 3)
	{
		Ray rRay = Ray(P, rVec, ray.r_IOR, bounces+1, giBounces, IS_REFLECT_RAY);
		HitInfo newHit;
		newHit.t = MIRO_TMAX;
		if (scene.trace(newHit, rRay, 0.001))
		{
			Vector3 reflection = newHit.obj->m_material->shade(rRay, newHit, scene);				
			Lr                += m_ks*m_reflectAmt*Rs*reflection;
			doEnv              = false;
		}
	}
	if (m_reflectAmt*Rs > 0.0f && doEnv)
	{
		Lr += m_ks * m_reflectAmt * Rs * getEnvironmentColor(rVec, scene);
	}

	if (m_refractAmt*Ts > 0.0f)
	{
		float snellsQ  = inIOR / outIOR;
		float sqrtPart = max(0.0f, sqrtf(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN)));
		Vector3 tVec   = (snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart)).normalized();	// Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf
		doEnv          = true;

		if (bounces < 3)																		// Do two bounces of refraction rays (2 bounces have already been used by reflection).
		{		
			ray.r_IOR.push(outIOR);
			Ray tRay = Ray(P, tVec, ray.r_IOR, bounces+1, giBounces, IS_REFRACT_RAY);			// Make a new refraction ray.
			HitInfo newHit;
			newHit.t = MIRO_TMAX;
			if (scene.trace(newHit, tRay, 0.001))
			{
				Vector3 refraction = newHit.obj->m_material->shade(tRay, newHit, scene);				
				Lt                += m_ks*m_refractAmt*Ts*refraction;
				doEnv              = false;
			}
		}
		if (doEnv)
		{
			Lt += m_ks * m_refractAmt * Ts * getEnvironmentColor(rVec, scene);
		}
	}

	Ld += m_ka;

	return (Ld + Ls) * (1.0f - Rs*m_reflectAmt - Ts*m_refractAmt) + Lr + Lt + m_Le; // Make sure we keep energy constant
}

/*Vector3 Blinn::shadeSSE(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
	union {float Ld[4];		__m128 _Ld; };
	union {float Lr[4];		__m128 _Lr; };
	union {float Lt[4];		__m128 _Lt; };
	union {float Ls[4];		__m128 _Ls; };
	union {float rayD[4];	__m128 _rayD; };
	union {float hitN[4];	__m128 _hitN; };
	union {float geoN[4];	__m128 _geoN; };
	union {float P[4];		__m128 _P; };
	union {float hitT[4];	__m128 _hitT; };
	union {};
	float u, v;

	const __m128 _zero			= setZero;
	const __m128 _neg			= setSSE(-1);
	const __m128 _neg0			= loadps(SSE_invertVec0);
	ALIGN_SSE Vector3 tmp1		= m_ks;
	ALIGN_SSE Vector3 tmp2		= m_kd;
	ALIGN_SSE Vector3 tmp3		= m_kt;
	ALIGN_SSE Vector3 tmp4		= m_ka;
	const __m128 _m_ks			= loadps(&tmp1.x);
	const __m128 _m_kd			= loadps(&tmp2.x);
	const __m128 _m_kt			= loadps(&tmp3.x);
	const __m128 _m_ka			= loadps(&tmp4.x);
	const __m128 _m_ior			= setSSE(m_ior);
	const __m128 _m_specExp		= setSSE(m_specExp);
	const __m128 _m_specAmt		= setSSE(m_specAmt);
	const __m128 _m_reflectAmt	= setSSE(m_reflectAmt);
	const __m128 _m_refractAmt	= setSSE(m_refractAmt);
	_Ld = setZero;
	_Lr = setZero;
	_Lt = setZero;
	_Ls = setZero;
	_hitT = setSSE(hit.t);

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = hit.tri->m_mesh;
	u_int meshIndex = hit.tri->m_index;

	ti3 = theMesh->m_vertexIndices[meshIndex];				// Get the geometric normal
	Vector3 geoNVec = cross(theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x], theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x]);
	for (int i = 0; i < 3; i++)
	{
		geoN[i] = geoNVec[i];
	}
	geoN[3] = 0.;
	_geoN = mulps(_geoN, fastrsqrtps(dotps(_geoN, _geoN, 0xFF)));

	ti3 = theMesh->m_normalIndices[meshIndex];				// Get the interpolated normal
	float a = hit.a; float b = hit.b; float c = 1.0f-a-b; 
	Vector3 nVec = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b));
	for (int i = 0; i < 3; i++)
	{
		hitN[i] = nVec[i];
	}
	hitN[3] = 0.;
	_hitN = mulps(_hitN, fastrsqrtps(dotps(_hitN, _hitN, 0xFF)));

	if (theMesh->m_texCoordIndices)							// If possible, get the interpolated u, v coordinates
	{
		ti3 = theMesh->m_texCoordIndices[meshIndex];
		u = theMesh->m_texCoords[ti3.x].x*c+theMesh->m_texCoords[ti3.y].x*a+theMesh->m_texCoords[ti3.z].x*b;
		v = theMesh->m_texCoords[ti3.x].y*c+theMesh->m_texCoords[ti3.y].y*a+theMesh->m_texCoords[ti3.z].y*b;
	}
	else
	{
		u = a;
		v = b;
	}

	_rayD = loadps(&ray.dx);
	const __m128 _viewDir = mulps(_rayD, _neg0);

	_P = addps(loadps(&ray.ox), mulps(_hitT, _rayD));		// Hit position

	bool flip = false;
	__m128 _vDotN		= dotps(_viewDir, _hitN, 0xFF);						// Find if the interpolated normal points back while
	__m128 _vDotGeoN	= dotps(_viewDir, _geoN, 0xFF);						// geometric normal is fine.
	bool nEqGeoN		= movemaskps(_vDotN) == movemaskps(_vDotGeoN);		// True if both signs are equal
	__m128 _normal		= nEqGeoN ? _hitN : _geoN;							// Use geometric normal if interpolated one points away
	_vDotN				= nEqGeoN ? _vDotN : _vDotGeoN;						// from object

	if (cmplessss_i(_vDotN, _zero))
	{
		flip = true;
		_vDotN = mulps(_vDotN, _neg);
		_normal = mulps(_normal, _neg);
	}

	const __m128 _rVec = addps(_rayD, mulps(setSSE(2), mulps(_vDotN, _normal)));			// get the reflection vector

	const Lights *lightlist = scene.lights();

	// loop over all of the lights
	Lights::const_iterator lightIter;
	for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
	{
		PointLight* pLight = *lightIter;

		union {float lDir[4]; __m128 _lDir; };
		ALIGN_SSE Vector3 lPos		= pLight->position();
		ALIGN_SSE Vector3 lColor	= pLight->color();
		ALIGN_SSE float	lWattage	= pLight->wattage();
		const __m128 _lPos			= loadps(&lPos.x);
		const __m128 _lColor		= loadps(&lColor.x);
		const __m128 _lWattage		= setSSE(lWattage);
		_lDir						= subps(_lPos, _P);

		// the inverse-squared falloff
		union {float distance[4]; __m128 _distance; };
		__m128 _falloff				= dotps(_lDir, _lDir, 0xFF);		// Distance to light, squared
		__m128 _distanceReciprocal	= fastrsqrtps(_falloff);			// Reciprocal of the distance
		_distance					= recipps(_distanceReciprocal);		// The distance
		_falloff					= recipps(_falloff);				// Reciprocal of the falloff

		// normalize the light direction
		_lDir = mulps(_lDir, _distanceReciprocal);

		__m128 _E;
		__m128 _nDotL = dotps(_normal, _lDir, 0xFF);
		if (!movemaskps(_nDotL))							// Check to see greater than 0 (we check sign bit)
		{
			_E = mulps(_nDotL, mulps(_lColor, mulps(_lWattage, mulps(setSSE(0.25), mulps(_falloff, setSSE(piRecip))))));		// Light irradiance
			Ray sRay = Ray(Vector3(P[0],P[1],P[2]), Vector3(lDir[0], lDir[1], lDir[2]), 1.0f, 0, IS_SHADOW_RAY);								// Create shadow ray
			HitInfo newHit;
			newHit.t = distance[0];

			if (!scene.trace(newHit, sRay, 0.001))								// Check for shadows
			{
				ALIGN_SSE float rDotL;
				storess(dotps(_rVec, _lDir, 0xFF), &rDotL);
				ALIGN_SSE float specPow = pow(rDotL, m_specExp);
				_Ls = addps(_Ls, mulps(_m_ks, mulps(_m_specAmt, setSSE(specPow))));			// Calculate specular component
				_Ld = addps(_Ld, mulps(_E, addps(_m_kd, _Ls)));								// Total = Light*(Diffuse+Spec)
			}
		}						
	}
	_Ld = addps(_Ld, _m_ka);
	return Vector3(Ld[0], Ld[1], Ld[2]);
}*/