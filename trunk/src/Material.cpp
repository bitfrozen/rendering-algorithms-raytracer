#include "Material.h"
#include "Scene.h"

Material::Material() : m_envExposure(1.0f), m_envMap(NULL), m_colorMap(NULL), 
	m_normalMap(NULL), m_reflectMap(NULL), m_refractMap(NULL), m_specularMap(NULL),
	m_alphaMap(NULL), m_sampleEnv(true), m_translucency(0.0f), m_disperse(false)
{
}

Material::~Material()
{
}

void Material::getCosineDistributedSamples(const int threadID, const Vector3 &N, Vector3 &out)
{
	ALIGN_SSE float sqrte2recip;
	ALIGN_SSE float sqrte2;
	ALIGN_SSE float sqrt1_e2recip;
	ALIGN_SSE float sqrt1_e2;

	//make random ray based on the cosine distribution
	const float e1 = Scene::getRand(threadID);
	float e2 = Scene::getRand(threadID);
	e2 = (e2 > 0.99) ? 0.99 : e2;

	Vector3 u = cross( ((fabs(N.x) > 0.1) ? Vector3(0,1,0) : Vector3(1,0,0)), N).normalize();
	Vector3 v = cross( N, u);

	float _2_PI_e1 = 2*PI*e1;

#ifndef NO_SSE
	fastrsqrtss(setSSE(e2), sqrte2recip);
	fastrsqrtss(setSSE(fabsf(1.0f - e2)), sqrt1_e2recip);
	recipss(setSSE(sqrte2recip), sqrte2);
	recipss(setSSE(sqrt1_e2recip), sqrt1_e2);
#else
	sqrte2   = sqrtf(e2);
	sqrt1_e2 = sqrtf(1.0f - e2);
#endif

	out = ((cos(_2_PI_e1)*sqrte2)*u + (sin(_2_PI_e1)*sqrte2)*v + sqrt1_e2*N).normalize();
}

const Vector3 Material::getEnvironmentColor(const Vector3& direction, const Scene& scene) const
{
	Vector3 envColor;
	if (m_envMap != NULL || scene.getEnvMap() != NULL) 
	{
		//environment map lookup
		if (m_envMap != NULL)
		{
			envColor = m_envMap->getLookupXYZ3(direction) * m_envExposure;
		} 
		else
		{
			envColor = scene.getEnvMap()->getLookupXYZ3(direction) * scene.getEnvExposure();
		}
	}
	else
	{
		envColor = scene.getBGColor();
	}
	return envColor;
}