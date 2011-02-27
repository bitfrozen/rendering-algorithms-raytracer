#include "Lambert.h"
#include "Ray.h"
#include "Scene.h"
#include "SSE.h"
#include "TriangleMesh.h"
#include "PointLight.h"

using namespace std;

Lambert::Lambert(const Vector3 & kd, const Vector3 & ka) :
    m_kd(kd), m_ka(ka)
{
}

Lambert::~Lambert()
{
}

const Vector3 Lambert::shade(const unsigned int threadID, const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
    Vector3 L		= Vector3(0.0f, 0.0f, 0.0f);
	Vector3 rayD	= Vector3(ray.d[0],ray.d[1],ray.d[2]);		// Ray direction
	Vector3 viewDir	= -rayD;									// View direction
	float u, v;
	Vector3 N, geoN;
	Vector3 diffuseColor = m_kd;

	Vector3 P = ray.getPoint(hit.t);

	hit.getAllInfos(N, geoN, u, v);

	if (m_texture != NULL) 
	{
		Vector4 texCol = m_texture->getLookup(u, v);
		diffuseColor = Vector3(texCol.x, texCol.y, texCol.z);
	}

    const Lights *lightlist = scene.lights();
    
    // loop over all of the lights
    Lights::const_iterator lightIter;
    for (lightIter = lightlist->begin(); lightIter != lightlist->end(); lightIter++)
    {
		float discard;
		Vector3 lightPower = (*lightIter)->sampleLight(threadID, P, N, scene, 0, discard);		
		L += lightPower * diffuseColor;								// Calculate Diffuse component
    }
    
    // add the ambient component
    L += m_ka;
    
    return L;
}
