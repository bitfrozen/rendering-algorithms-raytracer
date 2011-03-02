#include "ProxyObject.h"
#include "Ray.h"
#include "TriangleMesh.h"

ProxyObject::ProxyObject(Objects* m, BVH* b, Matrix4x4& t)
{
	m_objects = m;
	m_Matrix = (ProxyMatrix*)_aligned_malloc(sizeof(ProxyMatrix), 16);
	*m_Matrix = ProxyMatrix(t);
	m_BVH = b;
	m_objectType = PROXY_OBJECT;
}

void ProxyObject::getAABB(AABB* outBox)
{
	AABB tmpBox = m_BVH->getAABB();
	tmpBox.bbMin = m_Matrix->m_transform.multiplyAndDivideByW(tmpBox.bbMin);
	tmpBox.bbMax = m_Matrix->m_transform.multiplyAndDivideByW(tmpBox.bbMax);
	*outBox = tmpBox;
}

AABB ProxyObject::getAABB()
{
	AABB tmpBox = m_BVH->getAABB();
	tmpBox.bbMin = m_Matrix->m_transform.multiplyAndDivideByW(tmpBox.bbMin);
	tmpBox.bbMax = m_Matrix->m_transform.multiplyAndDivideByW(tmpBox.bbMax);
	return tmpBox;
}

const bool ProxyObject::intersect(const unsigned int threadID, HitInfo &result, const Ray& ray, const float tMin)
{
	Vector3 newRayO = m_Matrix->m_inverse.multiplyAndDivideByW(ray._o);
	Vector3 newRayD = m_Matrix->m_inverse * ray._d;
	HitInfo newHit;
	newHit.t = result.t;
	Ray newRay = Ray(threadID, newRayO, newRayD, ray.time, ray.r_IOR, 0, 0, ray.bounces_flags);

	bool hit = m_BVH->intersect(threadID, newHit, newRay, tMin);

	if (hit)
	{
		result.a = newHit.a;
		result.b = newHit.b;
		result.t = newHit.t;
		result.obj = newHit.obj;
		result.m_proxy = this;
	}
	return hit;
}

void ProxyObject::renderGL()
{
	for (int i = 0; i < m_objects->size(); i+=10000)
	{
		Vector3 v0, v1, v2;
		TriangleMesh::TupleI3 ti3 = (*m_objects)[i]->m_mesh->m_vertexIndices[(*m_objects)[i]->m_index];// [m_index];
		v0 = m_Matrix->m_transform.multiplyAndDivideByW( (*m_objects)[i]->m_mesh->m_vertices[ti3.x] ); //vertex a of triangle
		v1 = m_Matrix->m_transform.multiplyAndDivideByW( (*m_objects)[i]->m_mesh->m_vertices[ti3.y] ); //vertex b of triangle
		v2 = m_Matrix->m_transform.multiplyAndDivideByW( (*m_objects)[i]->m_mesh->m_vertices[ti3.z] ); //vertex c of triangle

		glBegin(GL_TRIANGLES);
		glVertex3f(v0.x, v0.y, v0.z);
		glVertex3f(v1.x, v1.y, v1.z);
		glVertex3f(v2.x, v2.y, v2.z);
		glEnd();
	}	
}

void ProxyObject::preCalc()
{
	if ((*m_objects)[0]->m_mesh->doPreCalc)
	{
		Objects::const_iterator objectIter;
		for (objectIter = m_objects->begin(); objectIter != m_objects->end(); objectIter++)
		{
			(*objectIter)->preCalc();
		}
	}
}

void ProxyObject::setupProxy(TriangleMesh* mesh, Material* mat, Objects* m, BVH* b)
{
	int numObjs = mesh->m_numTris;
	Object* t = new Object[numObjs];

	int i = numObjs-1;
	while (i >= 0)
	{
		t[i].setMesh(mesh);
		t[i].setIndex(i);
		t[i].setMaterial(mat);
		m->push_back(&t[i]);
		i--;
	}

	b->build(m);
}