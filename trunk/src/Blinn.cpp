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
			 float reflectAmt, float refractAmt, float specGloss) :
    m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt), m_ior(ior), 
	m_specExp(specExp), m_specAmt(specAmt),
	m_reflectAmt(reflectAmt), m_refractAmt(refractAmt), m_specGloss(specGloss)
{	
	m_lightEmitted = 0.0f;
	m_Le = Vector3(0.0f);
}

Blinn::~Blinn()
{
}

Vector3 Blinn::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
/*#ifndef NO_SSE
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

#else*/
	Vector3 Ld		= Vector3(    0.0f,    0.0f,    0.0f);		// Diffuse	
	Vector3 Lr		= Vector3(    0.0f,    0.0f,    0.0f);		// Reflection
	Vector3 Lt		= Vector3(    0.0f,    0.0f,    0.0f);		// Transmission
	Vector3 Ls		= Vector3(    0.0f,    0.0f,    0.0f);		// Specular
	Vector3 rayD	= Vector3(ray.d[0],ray.d[1],ray.d[2]);		// Ray direction
	Vector3 viewDir	= -rayD;									// View direction
	float u, v;

	Vector3 P;
#ifndef NO_SSE
	storeps(addps(ray._o, mulps(setSSE(hit.t), ray._d)), &P.x);	// Hit position
#else
	P = Vector3(ray.o[0] + hit.t*ray.d[0], ray.o[1] + hit.t*ray.d[1], ray.o[2] + hit.t*ray.d[2]);
#endif

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = hit.obj->m_mesh;
	u_int meshIndex       = hit.obj->m_index;

	ti3           = theMesh->m_vertexIndices[meshIndex];						// Get the geometric normal
	Vector3 edge0 = theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x];
	Vector3 edge1 = theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x];
	Vector3 geoN  = cross(edge0, edge1).normalized();

	ti3       = theMesh->m_normalIndices[meshIndex];				// Get the interpolated normal
	float c   = 1.0f-hit.a-hit.b;
	float a   = hit.a; float b = hit.b;

	Vector3 N = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b).normalized());

	if (theMesh->m_texCoordIndices)							// If possible, get the interpolated u, v coordinates
	{
		ti3 = theMesh->m_texCoordIndices[meshIndex];
		u   = theMesh->m_texCoords[ti3.x].x*c+theMesh->m_texCoords[ti3.y].x*a+theMesh->m_texCoords[ti3.z].x*b;
		v   = theMesh->m_texCoords[ti3.x].y*c+theMesh->m_texCoords[ti3.y].y*a+theMesh->m_texCoords[ti3.z].y*b;
	}
	else
	{
		u = a;
		v = b;
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

	Vector3 rVec = (rayD + 2*(vDotN)*theNormal);									// get the reflection vector
	
	////////////////////PATH TRACE
	
	Vector3 randD, ru, rv, envColor;
	float _2PIe1;
	ALIGN_SSE float sqrte2recip;
	ALIGN_SSE float sqrte2;
	ALIGN_SSE float sqrt1_e2recip;
	ALIGN_SSE float sqrt1_e2;
	ALIGN_SSE float e1;
	ALIGN_SSE float e2;
	Ray randRay;
	HitInfo newHit;
	float theta, phi, LookUpU, LookUpV;
	if (scene.m_pathTrace) {
		//compute diffuse component using path tracing.
		if (m_lightEmitted > 0.0f) {
			Ld += m_kd*m_lightEmitted*m_Le;
		}
		int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);
		if (bounces < scene.m_maxBounces) {
			//make random ray based on the cosine distribution
			if (Scene::randsIdx >= 990000)
			{
				Scene::randsIdx = 0;
				Scene::genRands();
			}
			e1 = Scene::rands[Scene::randsIdx++];
			e2 = Scene::rands[Scene::randsIdx++];
			e2 = (e2 > 0.99) ? 0.99 : e2;
			ru = cross(theNormal,edge0).normalize();
			rv = cross(theNormal,ru);
			_2PIe1 = 2*PI*e1;
			
#ifndef NO_SSE
			fastrsqrtss(setSSE(e2), sqrte2recip);
			fastrsqrtss(setSSE(fabsf(1-e2)), sqrt1_e2recip);
			recipss(setSSE(sqrte2recip), sqrte2);
			recipss(setSSE(sqrt1_e2recip), sqrt1_e2);
#else
			distance      = sqrtf(falloff);
			distanceRecip = 1.0f / distance;
			falloff       = 1.0f / falloff;
#endif
			randD = (cos(_2PIe1)*sqrte2*ru + sin(_2PIe1)*sqrte2*rv + sqrt1_e2*N).normalize();
			randRay.set(P, randD, 1.0f, bounces + 1, IS_PRIMARY_RAY);
			newHit.t = MIRO_TMAX;
			//trace the random ray
			
			if (scene.trace(newHit, randRay, 0.001))								// Check for shadows
			{
				Ld += m_kd*newHit.obj->m_material->shade(randRay,newHit,scene);
			}
			else
			{
				if (m_envMap != NULL || g_scene->getEnvMap() != NULL) 
				{
					//environment map lookup
					theta   = atan2(randD.z, randD.x) + PI;
					phi     = acos(randD.y);
					LookUpU = theta * 0.5 * piRecip;
					LookUpV = 1.0 - (phi * piRecip);
					if (m_envMap != NULL)
					{
						envColor        = m_envMap->getLookup3(LookUpU, LookUpV)*m_envExposure;
						Ld              += m_kd*envColor;
					} 
					else
					{
						envColor        = g_scene->getEnvMap()->getLookup3(LookUpU, LookUpV)*g_scene->getEnvExposure();
						Ld              += m_kd*envColor;
					}
				}
				else
				{
					Ld += m_kd*g_scene->getBGColor();
				}
			}

			{

				// Directly sample area lights
				const Lights *arealightlist = scene.lights();
				// loop over all of the area lights
				Lights::const_iterator arealightIter;
				for (arealightIter = arealightlist->begin(); arealightIter != arealightlist->end(); arealightIter++)
				{
					RectangleLight* rLight = dynamic_cast<RectangleLight*>(*arealightIter);
					if (rLight == 0)
						continue;

					Vector3 D = rLight->getPointOnLight() - P;

					// Directly sample square light, this would be rolled into the next loop over all lights...
					// I use the same two random numbers again, it shouldn't be an issue...
					//Vector3 p1, p2, p3;
					//p1 = Vector3(4.0, 5.5, -1.5);
					//p2 = Vector3(4.0, 5.5, -4.0);
					//p3 = Vector3(1.5, 5.5, -1.5);
					//Vector3 D = p1-P + e1*(p2-p1) + e2*(p3-p1); // Linear interpolation over the square
					
					// the inverse-squared falloff
					ALIGN_SSE float falloff = D.length2();
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
					D *= distanceRecip;

					float nDotL = max(0.0f, dot(theNormal, D));
					Vector3 E   = (nDotL * rLight->wattage()) * (0.25*falloff*piRecip);				// Light irradiance
					Ray sRay    = Ray(P, D, 1.0f, 0, IS_SHADOW_RAY);					// Create shadow ray

					if (nDotL)
					{
						HitInfo newHit;
						newHit.t = distance;
						if (scene.trace(newHit, sRay, 0.001))							// Check for shadows
						{
							E = 0;
						}
					}
					Ld += E*(m_kd + Ls);
				}
			}		

			// Directly sample pointlights
			const Lights *lightlist = scene.lights();

			// loop over all of the lights
			Lights::const_iterator lightIter;
			for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
			{
				PointLight* pLight = dynamic_cast<PointLight*>(*lightIter);
				if (pLight == 0)
					continue;

				Vector3 l = pLight->position() - P;

				// the inverse-squared falloff
				ALIGN_SSE float falloff = l.length2();
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

				// normalize the light direction
				l *= distanceRecip;

				float nDotL = max(0.0f, dot(theNormal, l));
				Vector3 E   = (nDotL * pLight->color() * pLight->wattage()) * (0.25*falloff*piRecip);		// Light irradiance
				Ray sRay    = Ray(P, l, 1.0f, 0, IS_SHADOW_RAY);											// Create shadow ray

				if (nDotL)
				{
					HitInfo newHit;
					newHit.t = distance;
					if (scene.trace(newHit, sRay, 0.001))								// Check for shadows
					{
						E = 0;
					}
				}

				Ls += m_ks*m_specAmt*pow(dot(rVec, l), m_specExp);								// Calculate specular component.
				Ld += E*(m_kd + Ls);															// Total = Light*(Diffuse+Spec)
			}

		} else {
			return Ld; //return no color since we did not hit a light
		}
 //////////////////////NO PATH TRACE
	} else { //no path tracing 
		const Lights *lightlist = scene.lights();

		// loop over all of the lights
		Lights::const_iterator lightIter;
		for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
		{
			PointLight* pLight = dynamic_cast<PointLight*>(*lightIter);
			if (pLight == 0)
				continue;

			Vector3 l = pLight->position() - P;

			// the inverse-squared falloff
			ALIGN_SSE float falloff = l.length2();
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

			// normalize the light direction
			l *= distanceRecip;

			float nDotL = max(0.0f, dot(theNormal, l));
			Vector3 E   = (nDotL * pLight->color() * pLight->wattage()) * (0.25*falloff*piRecip);		// Light irradiance
			Ray sRay    = Ray(P, l, 1.0f, 0, IS_SHADOW_RAY);											// Create shadow ray

			if (nDotL)
			{
				HitInfo newHit;
				newHit.t = distance;
				if (scene.trace(newHit, sRay, 0.001))								// Check for shadows
				{
					E = 0;
				}
			}

			Ls += m_ks*m_specAmt*pow(dot(rVec, l), m_specExp);								// Calculate specular component
			Ld += E*(m_kd + Ls);															// Total = Light*(Diffuse+Spec)
		}
	}  ///////////END

	// Here we do glossy reflection, however it de-coheres the rays too quickly
	if (m_specGloss < 1.0)
	{
		rVec = (m_specGloss*rVec + (1-m_specGloss)*randD).normalized();						// get the (randomized) reflection vector
	}

	float outIOR;									// outIOR is the IOR of the medium the ray will go into.
	float inIOR = ray.r_IOR();						// inIOR is the IOR of the medium the ray is currently in.
	if (flip)										// If this is true, then we've hit a back-facing polygon:
	{												// we're going out of the current material.
		ray.r_IOR.pop();							// Pop the back of the IORHistory so we get the IOR right, be careful and make sure this is OK.
		outIOR = ray.r_IOR();
	}
	else											// Normal. We're going (if possible) into a new material,
	{												// record its IOR in the ray's IORHistory
		outIOR = m_ior;
	}
	float Rs = 0; float Ts = 0;
	if (m_reflectAmt > 0.0 || m_refractAmt > 0.0)
	{
		Rs = fresnel(inIOR, outIOR, vDotN);												// Compute Fresnel's equations http://www.bramz.net/data/writings/reflection_transmission.pdf
		Ts = 1.0f-Rs;																	// Energy conservation
	}
	
	int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);						// Find out how many bounces we've taken
	bool doEnv  = true;

	if (m_reflectAmt*Rs > 0.0f && bounces < 4)											// Do two bounces of reflection rays.
	{
		Ray rRay = Ray(P, rVec, ray.r_IOR, bounces+1, IS_REFLECT_RAY);
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
		if (m_envMap != NULL || g_scene->getEnvMap() != NULL) 
		{
			//environment map lookup
			theta   = atan2(rVec.z, rVec.x) + PI;
			phi     = acos(rVec.y);
			LookUpU = theta * 0.5 * piRecip;
			LookUpV = 1.0 - (phi * piRecip);
			if (m_envMap != NULL)
			{
				Vector3 envColor = m_envMap->getLookup3(LookUpU, LookUpV);
				envColor        *= m_envExposure;
				Lr              += m_ks*m_reflectAmt*Rs*envColor;
			} 
			else
			{
				Vector3 envColor = g_scene->getEnvMap()->getLookup3(LookUpU, LookUpV);
				envColor        *= g_scene->getEnvExposure();
				Lr              += m_ks*m_reflectAmt*Rs*envColor;
			}
		}
		else
		{
			Lr += m_ks*m_reflectAmt*Rs*g_scene->getBGColor();
		}
	}

	if (m_refractAmt*Ts > 0.0f)
	{
		float snellsQ  = inIOR / outIOR;
		float sqrtPart = max(0.0f, sqrtf(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN)));
		Vector3 tVec   = (snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart)).normalized();	// Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf
		doEnv          = true;

		if (bounces < 4)																		// Do two bounces of refraction rays (2 bounces have already been used by reflection).
		{		
			ray.r_IOR.push(outIOR);
			Ray tRay = Ray(P, tVec, ray.r_IOR, bounces+1, IS_REFRACT_RAY);						// Make a new refraction ray.
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
			if (m_envMap != NULL || g_scene->getEnvMap() != NULL) 
			{
				//environment map lookup
				theta   = atan2(tVec.z, tVec.x) + PI;
				phi     = acos(tVec.y);
				LookUpU = theta * 0.5 * piRecip;
				LookUpV = 1.0 - (phi * piRecip);
				if (m_envMap != NULL)
				{
					Vector3 envColor = m_envMap->getLookup3(LookUpU, LookUpV);
					envColor        *= m_envExposure;
					Lt              += m_ks*m_refractAmt*Ts*envColor;
				} 
				else
				{
					Vector3 envColor = g_scene->getEnvMap()->getLookup3(LookUpU,LookUpV);
					envColor        *= g_scene->getEnvExposure();
					Lt              += m_ks*m_refractAmt*Ts*envColor;
				}
			}
			else
			{
				Lt += m_ks*m_refractAmt*Ts*g_scene->getBGColor();
			}
		}
	}

	Ld += m_ka;
	return Ld*(1.0-Rs*m_reflectAmt-Ts*m_refractAmt)+Lr+Lt;
//#endif
}

// Very quick and dirty image based lighting... still kinda cool!

/*

Vector3 Blinn::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
	Vector3 Ld		= Vector3(0.0f, 0.0f, 0.0f);			// Diffuse	
	Vector3 Lr		= Vector3(0.0f, 0.0f, 0.0f);			// Reflection
	Vector3 Lt		= Vector3(0.0f, 0.0f, 0.0f);			// Transmission
	Vector3 Ls		= Vector3(0.0f, 0.0f, 0.0f);			// Specular
	Vector3 rayD	= Vector3(ray.d[0], ray.d[1], ray.d[2]);		// Ray direction
	Vector3 viewDir	= -rayD;								// View direction
	float u, v;

	Vector3 P;
#ifndef NO_SSE
	storeps(addps(ray._o, mulps(setSSE(hit.t), ray._d)), &P.x);	// Hit position
#else
	P = Vector3(ray.ox + hit.t*ray.dx, ray.oy + hit.t*ray.dy, ray.oz + hit.t*ray.dz);
#endif

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = hit.obj->m_mesh;
	u_int meshIndex = hit.obj->m_index;

	ti3 = theMesh->m_vertexIndices[meshIndex];				// Get the geometric normal
	Vector3 edge0 = theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x];
	Vector3 edge1 = theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x];
	Vector3 geoN = cross(edge0, edge1).normalized();

	ti3 = theMesh->m_normalIndices[meshIndex];				// Get the interpolated normal
	float c = 1.0f-hit.a-hit.b;
	float a = hit.a; float b = hit.b;
	Vector3 N = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b).normalized());

	Vector3 T = cross(N, edge0).normalized();
	Vector3 biT = cross(T, N).normalized();

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

	float vDotN			= dot(viewDir, N);				// Find if the interpolated normal points back while
	float vDotGeoN		= dot(viewDir, geoN);			// geometric normal is fine.
	bool nEqGeoN		= (vDotN*vDotGeoN >= 0.0);		// True if both are equal
	Vector3 theNormal	= nEqGeoN ? N : geoN;			// Use geometric normal if interpolated one points away
	vDotN				= nEqGeoN ? vDotN : vDotGeoN;	// from object
	bool flip = false;
	if (vDotN < 0.0)
	{
		flip = true;
		vDotN = -vDotN;
		theNormal = -theNormal;
	}

	float num_samples_rec = 1.0f / num_samples;
	Ray newRay;
	HitInfo sHit;
	float theta, phi, lookupU, lookupV;

	for (int i = 0; i < num_samples; i++)
	{
		if (randsIdx >= 900000)
		{
			randsIdx = 0;
			genRands();
		}
		float sinTheta = sqrtf(rands[randsIdx++]);
		float newVecTheta = asinf(sinTheta);
		float newVecPhi   = 2*PI*rands[randsIdx++];
		//Vector3 newVec = sinTheta*cosf(newVecPhi)*biT + sinTheta*sinf(newVecPhi)*T + cosf(newVecTheta)*N;
		Vector3 newVec = Vector3(rands[randsIdx++], rands[randsIdx++], rands[randsIdx++]).normalized();
		//Vector3 newVec = Vector3(sinTheta*cosf(newVecPhi), cosf(newVecTheta), sinTheta*sinf(newVecPhi)).normalized();

		float E = 1;
		float newVdotN = dot(theNormal, newVec);
		newVec = (newVdotN > 0.0) ? newVec : -newVec;
		newVdotN = fabsf(newVdotN);

		newRay.set(P, newVec, 1.001f, 0, IS_SHADOW_RAY); 

		sHit.t = MIRO_TMAX;

		if (scene.trace(sHit, newRay, 0.001))								// Check for shadows
		{
			E = 0;
		}

		//environment map lookup
		theta = atan2(newVec.z, newVec.x) + PI;
		phi = acos(newVec.y);
		lookupU = theta * 0.5 * piRecip;
		lookupV = 1.0 - (phi * piRecip);
		Vector4 tmp = g_scene->getEnvMap()->getLookup3(lookupU, lookupV);
		Vector3 envColor = Vector3(tmp.x*num_samples_rec, tmp.y*num_samples_rec, tmp.z*num_samples_rec);
		Ld += m_kd*envColor*E*newVdotN;
	}
	
	Vector3 rVec = (rayD + 2*(vDotN)*theNormal);				// get the reflection vector

	/*const Lights *lightlist = scene.lights();

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

		float nDotL = max(0.0f, dot(theNormal, l));
		Vector3 E = (nDotL * pLight->color() * pLight->wattage()) / (4*falloff*PI);		// Light irradiance
		Ray sRay = Ray(P, l, 1.0f, 0, IS_SHADOW_RAY);									// Create shadow ray
		HitInfo sHit;
		sHit.t = MIRO_TMAX;
		if (nDotL)
		{
			if (scene.trace(sHit, sRay, 0.001))								// Check for shadows
			{
				E *= (sHit.t < distance) ? 0 : 1;
			}
		}

		Ls += m_ks*m_specAmt*pow(dot(rVec, l), m_specExp);								// Calculate specular component
		Ld += E*(m_kd + Ls);															// Total = Light*(Diffuse+Spec)
	}*/

	/*Ray::IORList IORHistory = ray.r_IOR;
	float outIOR;									// outIOR is the IOR of the medium the ray will go into.
	float inIOR = IORHistory.back();				// inIOR is the IOR of the medium the ray is currently in.
	if (flip)										// If this is true, then we've hit a back-facing polygon:
	{												// we're going out of the current material.
		if (IORHistory.size() > 1)					// Pop the back of the IORHistory so we get the IOR right, be careful and make sure this is OK.
		{
			IORHistory.pop_back();
		}		
		outIOR = IORHistory.back();
	}
	else											// Normal. We're going (if possible) into a new material,
	{												// record its IOR in the ray's IORHistory
		outIOR = m_ior;
	}
	float Rs = Material::fresnel(inIOR, outIOR, vDotN);									// Compute Fresnel's equations  http://www.bramz.net/data/writings/reflection_transmission.pdf
	float Ts = 1.0f-Rs;																	// Energy conservation
	int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);						// Find out how many bounces we've taken
	bool doEnv = true;

	if (m_reflectAmt*Rs > 0.0f && bounces < 3)											// Do two bounces of reflection rays.
	{
		Ray rRay = Ray(P, rVec, 
		ray.r_IOR, bounces+1, IS_REFLECT_RAY);
		HitInfo rHit;
		rHit.t = MIRO_TMAX;
		if (scene.trace(rHit, rRay, 0.001))
		{
			Vector3 reflection = rHit.obj->m_material->shade(rRay, rHit, scene);				
			Lr += m_ks*m_reflectAmt*Rs*reflection;
			doEnv = false;
		}
	}
	if (m_reflectAmt*Rs > 0.0f && doEnv)
	{
		if (m_envMap != NULL || g_scene->getEnvMap() != NULL) 
		{
			//environment map lookup
			theta = atan2(rVec.z, rVec.x) + PI;
			phi = acos(rVec.y);
			lookupU = theta * 0.5 * piRecip;
			lookupV = 1.0 - (phi * piRecip);
			if (m_envMap != NULL)
			{
				Vector3 envColor = m_envMap->getLookup3(lookupU, lookupV);
				envColor /= m_envExposure;
				Lr += m_ks*m_reflectAmt*Rs*envColor;
			} 
			else
			{
				Vector3 envColor = g_scene->getEnvMap()->getLookup3(lookupU, lookupV);
				envColor /= g_scene->getEnvExposure();
				Lr += m_ks*m_reflectAmt*Rs*envColor;
			}
		}
		else
		{
			Lr += m_ks*m_reflectAmt*Rs*g_scene->getBGColor();
		}
	}

	/*float snellsQ = inIOR / outIOR;
	float sqrtPart = max(0.0f, sqrt(1.0f - (snellsQ*snellsQ) * (1.0f-vDotN*vDotN)));
	Vector3 tVec = (snellsQ*rayD + theNormal*(snellsQ*vDotN - sqrtPart)).normalized();  // Get the refraction ray direction. http://www.bramz.net/data/writings/reflection_transmission.pdf
	doEnv = true;

	if (m_refractAmt*Ts > 0.0f && bounces < 3)											// Do two bounces of refraction rays (2 bounces have already been used by reflection).
	{		
		IORHistory.push_back(outIOR);
		Ray tRay = Ray(P, tVec, IORHistory, bounces+1, IS_REFRACT_RAY);					// Make a new refraction ray.
		HitInfo tHit;
		if (scene.trace(tHit, tRay, 0.001, MIRO_TMAX))
		{
			Vector3 refraction = tHit.material->shade(tRay, tHit, scene);				
			Lt += m_ks*m_refractAmt*Ts*refraction;
			doEnv = false;
		}
	}
	if (m_refractAmt*Ts > 0.0f && doEnv)
	{
		if (m_envMap != NULL || g_scene->getEnvMap() != NULL) 
		{
			//environment map lookup
			float theta = atan2(tVec.z, tVec.x) + PI;
			float phi = acos(tVec.y);
			float u = theta * 0.5 * piRecip;
			float v = 1.0 - (phi * piRecip);
			if (m_envMap != NULL)
			{
				Vector3 envColor = m_envMap->getLookup3(u,v);
				envColor /= m_envExposure;
				Lt += m_ks*m_refractAmt*Ts*envColor;
			} 
			else
			{
				Vector3 envColor = g_scene->getEnvMap()->getLookup3(u,v);
				envColor /= g_scene->getEnvExposure();
				Lt += m_ks*m_refractAmt*Ts*envColor;
			}
		}
		else
		{
			Lt += m_ks*m_refractAmt*Ts*g_scene->getBGColor();
		}
	}

	Ld += m_ka;*/

	/*return Ld*(1.0-Rs*m_reflectAmt-Ts*m_refractAmt)+Lr+Lt;
}*/