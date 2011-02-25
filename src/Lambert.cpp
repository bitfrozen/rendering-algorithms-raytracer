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

Vector3
Lambert::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
{
    Vector3 L = Vector3(0.0f, 0.0f, 0.0f);
	float u, v;
    
	Vector3 rayD = Vector3(ray.d[0], ray.d[1], ray.d[2]);
	Vector3 P;
#ifndef NO_SSE
	storeps(addps(ray._o, mulps(setSSE(hit.t), ray._d)), &P.x);	// Hit position
#else
	P = Vector3(ray.o[0] + hit.t*ray.d[0], ray.o[1] + hit.t*ray.d[1], ray.o[2] + hit.t*ray.d[2]);
#endif

	TriangleMesh::TupleI3 ti3;
	TriangleMesh* theMesh = hit.obj->m_mesh;
	u_int meshIndex = hit.obj->m_index;

	ti3 = theMesh->m_vertexIndices[meshIndex];
	Vector3 edge0 = theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x];
	Vector3 edge1 = theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x];
	Vector3 geoN = cross(edge0, edge1).normalized();

	ti3 = theMesh->m_normalIndices[meshIndex];
	float c = 1.0f-hit.a-hit.b;
	float a = hit.a; float b = hit.b;
	Vector3 N = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b).normalized());

	if (theMesh->m_texCoordIndices)
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
        float falloff = l.length2();
        
        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = max(0.0f, dot(N, l));
        Vector3 result = pLight->color();
        result *= m_kd;
        
        L += max(0.0f, nDotL/falloff * pLight->wattage() / PI) * result;
    }
    
	if (m_texture != NULL) {
		Vector4 texCol = m_texture->getLookup(u, v);
		L = L*Vector3(texCol.x, texCol.y, texCol.z);
	}
    // add the ambient component
    L += m_ka;
    
    return L;
}
