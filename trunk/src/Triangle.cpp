#include "Triangle.h"
#include "TriangleMesh.h"
#include "Ray.h"

Triangle::Triangle(TriangleMesh * m, unsigned int i) :
    m_mesh(m), m_index(i)
{

}


Triangle::~Triangle()
{

}

void
Triangle::preCalc()
{
	m_mesh->preCalc();
}

void
Triangle::renderGL()
{
	Vector3 v0, v1, v2;
	TriangleMesh::TupleI3 ti3;
    for (int i = 0; i < m_mesh->numTris(); i++)
	{
		ti3 = m_mesh->vIndices()[i];// [m_index];
		v0 = m_mesh->vertices()[ti3.x]; //vertex a of triangle
		v1 = m_mesh->vertices()[ti3.y]; //vertex b of triangle
		v2 = m_mesh->vertices()[ti3.z]; //vertex c of triangle

		glBegin(GL_TRIANGLES);
			glVertex3f(v0.x, v0.y, v0.z);
			glVertex3f(v1.x, v1.y, v1.z);
			glVertex3f(v2.x, v2.y, v2.z);
		glEnd();
	}
}

bool
Triangle::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
	if (m_mesh->intersect(result, r, tMin, tMax))
	{
		result.material = this->m_material;
		return true;
	}
	return false;
}