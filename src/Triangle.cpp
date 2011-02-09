#include "TriangleMesh.h"
#include "Triangle.h"
#include "Ray.h"

Triangle::Triangle(TriangleMesh* m, u_int i) :
    m_mesh(m), m_index(i)
{
	updateAABB();
}
	
Triangle::~Triangle()
{
}

void Triangle::updateAABB()
{
	if (m_mesh)
	{
		Vector3 A, B, C;
		TriangleMesh::TupleI3 ti3 = m_mesh->m_vertexIndices[m_index];// [m_index];
		A = m_mesh->m_vertices[ti3.x]; //vertex a of triangle
		B = m_mesh->m_vertices[ti3.y]; //vertex b of triangle
		C = m_mesh->m_vertices[ti3.z]; //vertex c of triangle
		*m_bBox = AABB(Vector3(std::min(A.x, std::min(B.x, C.x)),
			std::min(A.y, std::min(B.y, C.y)),
			std::min(A.z, std::min(B.z, C.z))),
			Vector3(std::max(A.x, std::max(B.x, C.x)),
			std::max(A.y, std::max(B.y, C.y)),
			std::max(A.z, std::max(B.z, C.z))));
	}
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
	TriangleMesh::TupleI3 ti3 = m_mesh->m_vertexIndices[m_index];// [m_index];
	v0 = m_mesh->m_vertices[ti3.x]; //vertex a of triangle
	v1 = m_mesh->m_vertices[ti3.y]; //vertex b of triangle
	v2 = m_mesh->m_vertices[ti3.z]; //vertex c of triangle
	
	glBegin(GL_TRIANGLES);
		glVertex3f(v0.x, v0.y, v0.z);
		glVertex3f(v1.x, v1.y, v1.z);
		glVertex3f(v2.x, v2.y, v2.z);
	glEnd();
}

bool
Triangle::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
	if (m_mesh->intersect(result, r, tMin, tMax, m_index))
	{
		result.material = this->m_material;
		return true;	
	}
	return false;
}