#include "Blinn.h"
#include "Ray.h"
#include "Scene.h"
#include <time.h>
#include <math.h>
#include "Miro.h"
#include "SSE.h"
#include "TriangleMesh.h"

using namespace std;

Blinn::Blinn(const Vector3 & kd, const Vector3 & ka,
			 const Vector3 & ks, const Vector3 & kt,
			 float ior, float specExp, float specAmt,
			 float reflectAmt, float refractAmt) :
    m_kd(kd), m_ka(ka), m_ks(ks), m_kt(kt), m_ior(ior), 
	m_specExp(specExp), m_specAmt(specAmt),
	m_reflectAmt(reflectAmt), m_refractAmt(refractAmt)
{	
	//genRands();
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
	
	if (scene.m_pathTrace) {
		//compute diffuse component using path tracing.
		if (m_lightEmitted > 0.0f) {
			//insert BRDF here if we plan on using one
			Ld += m_lightEmitted*m_Le;
			return Ld;
		}
		int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);
		if (bounces < scene.m_maxBounces) {
			//make random ray based on the cosine distribution
			MTRand drand;
			float e1 = drand();
			float e2 = drand();
			Vector3 rw = geoN;
			Vector3 ru = edge0;
			ru = cross(rw,ru);
			ru = ru.normalize();
			Vector3 rv = cross(rw,ru);
			Vector3 randD = cos(2*PI*e1)*sqrt(e2)*ru + sin(2*PI*e1)*sqrt(e2)*rv + sqrt(1 - e2)*rw;
			Ray randRay = Ray(P, randD, 1.0f, bounces + 1, IS_PRIMARY_RAY);

			//trace the random ray
			HitInfo newHit;
			if (scene.trace(newHit, randRay, 0.001))								// Check for shadows
			{
				Ld += m_kd*newHit.obj->m_material->shade(randRay,newHit,scene);
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
			PointLight* pLight = *lightIter;

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

	Ray::IORList IORHistory = ray.r_IOR;
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
	float Rs = 0; float Ts = 0;
	if (m_reflectAmt > 0.0 || m_refractAmt > 0.0)
	{
		Rs = fresnel(inIOR, outIOR, vDotN);												// Compute Fresnel's equations  http://www.bramz.net/data/writings/reflection_transmission.pdf
		Ts = 1.0f-Rs;																	// Energy conservation
	}
	
	int bounces = GET_BOUNCES(ray.bounces_flags & BOUNCES_MASK);						// Find out how many bounces we've taken
	bool doEnv  = true;

	if (m_reflectAmt*Rs > 0.0f && bounces < 10)											// Do two bounces of reflection rays.
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
			float theta   = atan2(rVec.z, rVec.x) + PI;
			float phi     = acos(rVec.y);
			float LookUpU = theta * 0.5 * piRecip;
			float LookUpV = 1.0 - (phi * piRecip);
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

		if (bounces < 10)																		// Do two bounces of refraction rays (2 bounces have already been used by reflection).
		{		
			IORHistory.push_back(outIOR);
			Ray tRay = Ray(P, tVec, IORHistory, bounces+1, IS_REFRACT_RAY);						// Make a new refraction ray.
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
				float theta   = atan2(tVec.z, tVec.x) + PI;
				float phi     = acos(tVec.y);
				float LookUpU = theta * 0.5 * piRecip;
				float LookUpV = 1.0 - (phi * piRecip);
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

/*int Blinn::randsIdx = 0;		// Current index into random number array
float Blinn::rands[1000000];	// Array of random numbers

void Blinn::genRands()			// Run the random number generator. This way we don't run into
{								// threading problems...
	MTRand_int32 drand(clock());
	for (int i = 0; i < 1000000; i++)
	{
		rands[i] = ((float)drand()+0.5) * IntRecip;
	}
}

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

	float vDotN			= dot(viewDir, N);					// Find if the interpolated normal points back while
	float vDotGeoN		= dot(viewDir, geoN);					// geometric normal is fine.
	bool nEqGeoN		= (vDotN*vDotGeoN >= 0.0);				// True if both are equal
	Vector3 theNormal	= nEqGeoN ? N : geoN;			// Use geometric normal if interpolated one points away
	vDotN				= nEqGeoN ? vDotN : vDotGeoN;			// from object
	bool flip = false;
	if (vDotN < 0.0)
	{
		flip = true;
		vDotN = -vDotN;
		theNormal = -theNormal;
	}

	float num_samples_rec = 1.0f / num_samples;
	Ray newRay;
	float theta, phi, lookupU, lookupV;

	for (int i = 0; i < num_samples; i++)
	{
		if (randsIdx >= 900000)
		{
			randsIdx = 0;
			genRands();
		}
		float newVecPhi = asinf(sqrtf(rands[randsIdx++]));
		float newVecTheta   = 2*PI*rands[randsIdx++];
		Vector3 newVec = Vector3(sinf(newVecTheta)*cosf(newVecPhi), sinf(newVecTheta)*sinf(newVecPhi), cosf(newVecTheta)).normalized();

		float E = 1;
		//float newVdotN = dot(theNormal, newVec);
		//newVec = (newVdotN > 0.0) ? newVec : -newVec;
		//newVdotN = fabsf(newVdotN);

		newRay.o[0] = P.x; newRay.o[1] = P.y; newRay.o[2] = P.z; newRay.o[3] = 1.0f;
		newRay.d[0] = newVec.x; newRay.d[1] = newVec.y; newRay.d[2] = newVec.z; newRay.d[3] = 0.0f;
		newRay.id[0] = 1.0f/newVec.x; newRay.id[1] = 1.0f/newVec.y; newRay.id[2] = 1.0f/newVec.z; newRay.id[3] = 0.0f;
		newRay.bounces_flags = IS_SHADOW_RAY;

		HitInfo sHit;
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
		Ld += m_kd*envColor*E;//*newVdotN;
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