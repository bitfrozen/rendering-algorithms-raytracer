#include "TriangleMesh.h"
#include "Triangle.h"
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

void Triangle::getBounds(Vector3& bMin, Vector3& bMax)
{
	getBounds(bMin, bMax, 0, m_mesh->numTris()-1);
}

void Triangle::getBounds(Vector3& bMin, Vector3& bMax, int iMin, int iMax)
{
	bMin = MIRO_TMAX;
	bMax = -MIRO_TMAX;
	Vector3 v0, v1, v2;
	TriangleMesh::TupleI3 ti3;
	for (int i = iMin; i <= iMax; i++)
	{
		ti3 = m_mesh->vIndices()[i];
		v0  = m_mesh->vertices()[ti3.x]; //vertex a of triangle
		v1  = m_mesh->vertices()[ti3.y]; //vertex b of triangle
		v2  = m_mesh->vertices()[ti3.z]; //vertex c of triangle
		bMin.x = std::min(bMin.x, std::min(v0.x, std::min(v1.x, v2.x)));
		bMin.y = std::min(bMin.y, std::min(v0.y, std::min(v1.y, v2.y)));
		bMin.z = std::min(bMin.z, std::min(v0.z, std::min(v1.z, v2.z)));
		bMax.x = std::max(bMax.x, std::max(v0.x, std::max(v1.x, v2.x)));
		bMax.y = std::max(bMax.y, std::max(v0.y, std::max(v1.y, v2.y)));
		bMax.z = std::max(bMax.z, std::max(v0.z, std::max(v1.z, v2.z)));
	}	
}