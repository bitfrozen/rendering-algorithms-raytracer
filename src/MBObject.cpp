#include "TriangleMesh.h"
#include "MBObject.h"
#include "Ray.h"
#include "SSE.h"
#include <omp.h>

MBObject::MBObject(Material* mat, TriangleMesh* m, TriangleMesh* m2, u_int i) : Object(mat, m, i)
{
	m_mesh_t2 = m2;
	m_objectType = MB_OBJECT;
}

MBObject::~MBObject()
{
}

void MBObject::preCalc()
{
	m_mesh->preCalc();
	m_mesh_t2->preCalc();
}

const bool MBObject::intersect(const unsigned int threadID, HitInfo& result, const Ray& r, const float tMin)
{
	float t = min(1.f, max(0.f, r.time));
	
	//#pragma omp atomic
	Ray::rayTriangleIntersections[threadID]++;

	bool hit = false;
	ALIGN_SSE float newT;
	bool shadow = r.bounces_flags & IS_SHADOW_RAY;
#ifndef NO_SSE													// Implementation of Havel-Herout(2010)

	// Get the triangle at time t
	const __m128 _time = setSSE(t);
	const __m128 _1_time = subps(_onesps, _time);
	const __m128 _n = addps(mulps(_time, m_mesh_t2->m_preCalcTris[m_index]._n), mulps(_1_time, m_mesh->m_preCalcTris[m_index]._n));
	const __m128 _u = addps(mulps(_time, m_mesh_t2->m_preCalcTris[m_index]._u), mulps(_1_time, m_mesh->m_preCalcTris[m_index]._u));
	const __m128 _v = addps(mulps(_time, m_mesh_t2->m_preCalcTris[m_index]._v), mulps(_1_time, m_mesh->m_preCalcTris[m_index]._v));

	const __m128 det = dotps(_n, r._d, 0x7f);
	const __m128 invertVec4 = loadps(SSE_invertVec);
	const __m128 dett = dotps(mulps(invertVec4, _n), r._o, 0xff);
	const __m128 oldt = loadss(&result.t);

	if ((movemaskps(xorps(dett,subss(mulss(oldt, det), dett)))&1) == 0)
	{
		const __m128 detp = addps(mulps(r._o, det), mulps(dett, r._d));
		const __m128 detu = dotps(detp, _u, 0xf1);
		if ((movemaskps(xorps(detu, subss(det, detu)))&1) == 0)
		{
			const __m128 detv = dotps(detp, _v, 0xf1);

			if ((movemaskps(xorps(detv, subss(det, addss(detu, detv))))&1) == 0)
			{
				__m128 inv_det = recipss(det);
				storess(mulss(dett, inv_det), &newT);
				if (newT >= tMin && newT < result.t)
				{
					result.t = newT;
					if (shadow) return true;

					storess(mulss(detu, inv_det), &result.a);
					storess(mulss(detv, inv_det), &result.b);
					result.obj = this;
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

	edge[0] = t*(m_mesh->m_vertices[ti3.y] - m_mesh->m_vertices[ti3.x]) + (1.f-t)*(m_mesh_t2->m_vertices[ti3.y] - m_mesh_t2->m_vertices[ti3.x]);
	edge[1] = t*(m_mesh->m_vertices[ti3.z] - m_mesh->m_vertices[ti3.x]) + (1.f-t)*(m_mesh_t2->m_vertices[ti3.z] - m_mesh_t2->m_vertices[ti3.x]);

	pvec = cross(rd, edge[1]);

	det = dot(edge[0], pvec);

#ifdef TEST_CULL
	tvec = ro - (t*m_mesh->m_vertices[ti3.x] + (1-t)*m_mesh_t2->m_vertices[ti3.x]);
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

	tvec = ro - (t*m_mesh->m_vertices[ti3.x] + (1-t)*m_mesh_t2->m_vertices[ti3.x]);
	a = dot(tvec, pvec) * inv_det;
	if (a < 0.0 || a > 1.0) return false;

	qvec = cross(tvec, edge[0]);
	b = dot(rd, qvec) * inv_det;
	if (b < 0.0 || a + b > 1.0) return false;

	newT = dot(edge[1], qvec)*inv_det;
#endif
	if (newT >= tMin && newT < result.t)
	{		
		result.t = newT;
		if (shadow) return true;
		result.a = a;
		result.b = b;
		result.obj = this;
		return true;
	}	
#endif
	return hit;
}

void MBObject::getAABB(AABB* outBox)
{
	*outBox = AABB(m_mesh->getAABB(m_index), m_mesh_t2->getAABB(m_index));
}

AABB MBObject::getAABB()
{
	return AABB(m_mesh->getAABB(m_index), m_mesh_t2->getAABB(m_index));
}