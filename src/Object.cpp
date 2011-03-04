#include "TriangleMesh.h"
#include "Object.h"
#include "Ray.h"
#include "SSE.h"
#include <omp.h>

Object::Object(Material* mat, TriangleMesh* m, u_int i) :
    m_material(mat), m_mesh(m), m_index(i)
{
	m_objectType = OBJECT;
}
	
Object::~Object()
{
}

void Object::preCalc()
{
	m_mesh->preCalc();
}

void Object::renderGL()
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

const bool Object::intersect(const unsigned int threadID, HitInfo &result, const Ray& r, const float tMin)
{
//#pragma omp atomic
	Ray::rayTriangleIntersections[threadID]++;

	bool hit = false;
	ALIGN_SSE float newT;
	bool shadow = r.bounces_flags & IS_SHADOW_RAY;
#ifndef NO_SSE													// Implementation of Havel-Herout(2010)

	const __m128 det = dotps(m_mesh->m_preCalcTris[m_index]._n, r._d, 0x7f);
	const __m128 invertVec4 = loadps(SSE_invertVec);
	const __m128 dett = dotps(mulps(invertVec4, m_mesh->m_preCalcTris[m_index]._n), r._o, 0xff);
	const __m128 oldt = loadss(&result.t);

	if ((movemaskps(xorps(dett,subss(mulss(oldt, det), dett)))&1) == 0)
	{
		const __m128 detp = addps(mulps(r._o, det), mulps(dett, r._d));
		const __m128 detu = dotps(detp, m_mesh->m_preCalcTris[m_index]._u, 0xf1);
		if ((movemaskps(xorps(detu, subss(det, detu)))&1) == 0)
		{
			const __m128 detv = dotps(detp, m_mesh->m_preCalcTris[m_index]._v, 0xf1);

			if ((movemaskps(xorps(detv, subss(det, addss(detu, detv))))&1) == 0)
			{
				__m128 inv_det = recipss(det);
				storess(mulss(dett, inv_det), &newT);
				if (newT >= tMin && newT < result.t)
				{
					if (m_material->m_alphaMap)
					{
						ALIGN_SSE float a;
						ALIGN_SSE float b;
						storess(mulss(detu, inv_det), &a);
						storess(mulss(detv, inv_det), &b);
						const float c = 1.0f-a-b;
						float uCoord, vCoord;
						if (m_mesh->m_texCoordIndices)			// If possible, get the interpolated u, v coordinates
						{
							const TriangleMesh::TupleI3 ti3 = m_mesh->m_texCoordIndices[m_index];
							uCoord = m_mesh->m_texCoords[ti3.x].x*c+m_mesh->m_texCoords[ti3.y].x*a+m_mesh->m_texCoords[ti3.z].x*b;
							vCoord = m_mesh->m_texCoords[ti3.x].y*c+m_mesh->m_texCoords[ti3.y].y*a+m_mesh->m_texCoords[ti3.z].y*b;
						}
						else													// We always return texture coordinates
						{
							uCoord = a;
							vCoord = b;
						}
						float alpha = m_material->m_alphaMap->getLookupAlpha(uCoord, vCoord);
						if (alpha < 0.5f) return false;
						result.t = newT;
						if (shadow) return true;

						result.a = a;
						result.b = b;
						result.obj = this;
						result.m_proxy = NULL;
						return true;
					}
					result.t = newT;
					if (shadow) return true;

					storess(mulss(detu, inv_det), &result.a);
					storess(mulss(detv, inv_det), &result.b);
					result.obj = this;
					result.m_proxy = NULL;
					return true;
				}
			}
		}
	}
#else
	Vector3 edge[2], tvec, pvec, qvec,rd,ro;				// Implementation of Moller-Trumbore
	rd = Vector3(r.d[0],r.d[1],r.d[2]);
	ro = Vector3(r.o[0],r.o[1],r.o[2]);
	float det, inv_det, a, b;
	TriangleMesh::TupleI3 ti3;
	ti3 = m_mesh->m_vertexIndices[m_index];
	edge[0] = m_mesh->m_vertices[ti3.y] - m_mesh->m_vertices[ti3.x];
	edge[1] = m_mesh->m_vertices[ti3.z] - m_mesh->m_vertices[ti3.x];

	pvec = cross(rd, edge[1]);

	det = dot(edge[0], pvec);

#ifdef TEST_CULL
	tvec = ro - m_mesh->m_vertices[ti3.x];
	a = dot(tvec, pvec);	
	if (a < 0.0 || a > det) return false;

	qvec = cross(tvec, edge[0]);
	b = dot(rd, qvec);
	if (b < 0.0 || a + b > det) return false;

	inv_det = 1.0 / det;
	float newT = dot(edge[1], qvec)*inv_det;
	a *= inv_det;
	b *= inv_det;
#else
	inv_det = 1.0 / det;

	tvec = ro - m_mesh->m_vertices[ti3.x];
	a = dot(tvec, pvec) * inv_det;
	if (a < 0.0 || a > 1.0) return false;

	qvec = cross(tvec, edge[0]);
	b = dot(rd, qvec) * inv_det;
	if (b < 0.0 || a + b > 1.0) return false;

	newT = dot(edge[1], qvec)*inv_det;
#endif
	if (newT >= tMin && newT < result.t)
	{	
		if (m_material->m_alphaMap)
		{
			const float c = 1.0f-a-b;
			float uCoord, vCoord;
			if (m_mesh->m_texCoordIndices)			// If possible, get the interpolated u, v coordinates
			{
				const TriangleMesh::TupleI3 ti3 = m_mesh->m_texCoordIndices[m_index];
				uCoord = m_mesh->m_texCoords[ti3.x].x*c+m_mesh->m_texCoords[ti3.y].x*a+m_mesh->m_texCoords[ti3.z].x*b;
				vCoord = m_mesh->m_texCoords[ti3.x].y*c+m_mesh->m_texCoords[ti3.y].y*a+m_mesh->m_texCoords[ti3.z].y*b;
			}
			else													// We always return texture coordinates
			{
				uCoord = a;
				vCoord = b;
			}
			float alpha = m_material->m_alphaMap->getLookupAlpha(uCoord, vCoord);
			if (alpha < 0.5f) return false;

			result.t = newT;
			if (shadow) return true;

			result.a = a;
			result.b = b;
			result.obj = this;
			result.m_proxy = NULL;
			return true;
		}
		result.t = newT;
		if (shadow) return true;
		result.a = a;
		result.b = b;
		result.obj = this;
		result.m_proxy = NULL;
		return true;
	}	
#endif
	return hit;
}

void Object::getAABB(AABB* outBox)
{
	m_mesh->getAABB(m_index, outBox);
}

AABB Object::getAABB()
{
	return m_mesh->getAABB(m_index);
}

int Object::sortByArea(const void* s1, const void* s2)
{
	float left, right;
	left = (*(Object**)s1)->getAABB().getArea();
	right = (*(Object**)s2)->getAABB().getArea();

	if (left < right)
	{
		return -1;
	}
	else if (left > right)
	{
		return 1;
	}
	return 0;
}

int Object::sortByXComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

	if (left.x < right.x)
	{
		return -1;
	}
	else if (left.x > right.x)
	{
		return 1;
	}
	return 0;
}

int Object::sortByYComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

	if (left.y < right.y)
	{
		return -1;
	}
	else if (left.y > right.y)
	{
		return 1;
	}
	return 0;
}

int Object::sortByZComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

	if (left.z < right.z)
	{
		return -1;
	}
	else if (left.z > right.z)
	{
		return 1;
	}
	return 0;
}