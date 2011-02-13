#include "TriangleMesh.h"
#include "Triangle.h"
#include "Scene.h"
#include <xmmintrin.h>
#include <smmintrin.h>

TriangleMesh::TriangleMesh() :
    m_normals(0),
    m_vertices(0),
    m_texCoords(0),
    m_normalIndices(0),
    m_vertexIndices(0),
    m_texCoordIndices(0),
	m_preCalcTris(0),
	doPreCalc(true)
{
}

TriangleMesh::~TriangleMesh()
{
    delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_normalIndices;
    delete [] m_vertexIndices;
    delete [] m_texCoordIndices;
	_aligned_free(m_preCalcTris);
}

void TriangleMesh::preCalc()
{
#ifndef NO_SSE
	if (doPreCalc)
	{
		if (m_preCalcTris) _aligned_free(m_preCalcTris);
		m_preCalcTris = (PrecomputedTriangle*)_aligned_malloc(sizeof(PrecomputedTriangle)*m_numTris, 16); //new PrecomputedTriangle[m_numTris];

		Vector3 A, B, C, N, N1, N2, AC, AB;
		float d, d1, d2;
		TriangleMesh::TupleI3 ti3;
		for (u_int i = 0; i < m_numTris; i++)
		{
			ti3 = m_vertexIndices[i];// [m_index];
			A = m_vertices[ti3.x]; //vertex a of triangle
			B = m_vertices[ti3.y]; //vertex b of triangle
			C = m_vertices[ti3.z]; //vertex c of triangle
			AC = C-A;
			AB = B-A;
			N = cross(AB,AC);
			d = dot(A,N);
			float Nsq = 1.0 / N.length2();

			N1 = cross(AC,N) * Nsq;
			d1 = -dot(N1,A);
			N2 = cross(N,AB) * Nsq;
			d2 = -dot(N2,A);		
			m_preCalcTris[i].nx = N.x;
			m_preCalcTris[i].ny = N.y;
			m_preCalcTris[i].nz = N.z;
			m_preCalcTris[i].nd = d;
			m_preCalcTris[i].ux = N1.x;
			m_preCalcTris[i].uy = N1.y;
			m_preCalcTris[i].uz = N1.z;
			m_preCalcTris[i].ud = d1;
			m_preCalcTris[i].vx = N2.x;
			m_preCalcTris[i].vy = N2.y;
			m_preCalcTris[i].vz = N2.z;
			m_preCalcTris[i].vd = d2;
		}
		doPreCalc = false;
	}
#endif
}

void TriangleMesh::getAABB(u_int index, AABB* outBox)
{
	Vector3 A, B, C, bbMin, bbMax;
	TriangleMesh::TupleI3 ti3;

	ti3 = m_vertexIndices[index];// [m_index];
	A = m_vertices[ti3.x]; //vertex a of triangle
	B = m_vertices[ti3.y]; //vertex b of triangle
	C = m_vertices[ti3.z]; //vertex c of triangle

	bbMin.x = std::min(A.x, std::min(B.x, C.x));
	bbMin.y = std::min(A.y, std::min(B.y, C.y));
	bbMin.z = std::min(A.z, std::min(B.z, C.z));
	bbMax.x = std::max(A.x, std::max(B.x, C.x));
	bbMax.y = std::max(A.y, std::max(B.y, C.y));
	bbMax.z = std::max(A.z, std::max(B.z, C.z));

	*outBox = AABB(bbMin, bbMax);
	outBox->doAC();
}

AABB TriangleMesh::getAABB(u_int index)
{
	Vector3 A, B, C, bbMin, bbMax;
	TriangleMesh::TupleI3 ti3;

	ti3 = m_vertexIndices[index];// [m_index];
	A = m_vertices[ti3.x]; //vertex a of triangle
	B = m_vertices[ti3.y]; //vertex b of triangle
	C = m_vertices[ti3.z]; //vertex c of triangle

	bbMin.x = std::min(A.x, std::min(B.x, C.x));
	bbMin.y = std::min(A.y, std::min(B.y, C.y));
	bbMin.z = std::min(A.z, std::min(B.z, C.z));
	bbMax.x = std::max(A.x, std::max(B.x, C.x));
	bbMax.y = std::max(A.y, std::max(B.y, C.y));
	bbMax.z = std::max(A.z, std::max(B.z, C.z));

	AABB outBox = AABB(bbMin, bbMax);
	outBox.doAC();
	return outBox;
}

bool TriangleMesh::intersect(HitInfo& result, const Ray& r, float tMin, float tMax, u_int index)
{
	bool hit = false;
	TriangleMesh::TupleI3 ti3;
	ALIGN_SSE float newT, a, b;
	bool shadow = r.bounces_flags & IS_SHADOW_RAY;
#ifndef NO_SSE													// Implementation of Havel-Herout(2010), works fine for diffuse shading but breaks up for reflection/refraction
	__declspec(align(16)) const float arr[4] = {-1,-1,-1,1};
	__declspec(align(16)) const float arr2[4] = {2,2,2,2};
	const __m128 int_coef = _mm_load_ps(arr);
	const __m128 coef_2 = _mm_load_ps(arr2);

	const __m128 o = _mm_load_ps(&r.ox);
	const __m128 d = _mm_load_ps(&r.dx);
	const __m128 n = _mm_load_ps(&m_preCalcTris[index].nx);

	const __m128 det = _mm_dp_ps(n, d, 0x7f);
	const __m128 dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
	const __m128 oldt = _mm_load_ss(&tMax);

	if ((_mm_movemask_ps(_mm_xor_ps(dett,_mm_sub_ss(_mm_mul_ss(oldt, det), dett)))&1) == 0)
	{
		const __m128 detp = _mm_add_ps(_mm_mul_ps(o, det), _mm_mul_ps(dett, d));
		const __m128 detu = _mm_dp_ps(detp, _mm_load_ps(&m_preCalcTris[index].ux), 0xf1);
		if ((_mm_movemask_ps(_mm_xor_ps(detu, _mm_sub_ss(det, detu)))&1) == 0)
		{
			const __m128 detv = _mm_dp_ps(detp, _mm_load_ps(&m_preCalcTris[index].vx), 0xf1);

			if ((_mm_movemask_ps(_mm_xor_ps(detv, _mm_sub_ss(det, _mm_add_ss(detu, detv))))&1) == 0)
			{
				__m128 inv_det = _mm_rcp_ss(det);
				inv_det = _mm_sub_ss(_mm_mul_ss(coef_2, inv_det), _mm_mul_ss(det, _mm_mul_ss(inv_det, inv_det)));		// One step of Newton-Raphson, otherwise accuracy is too low.
				_mm_store_ss(&newT, _mm_mul_ss(dett, inv_det));
				if (newT >= tMin && newT < result.t)
				{
					result.t = newT;
					if (!shadow)
					{
						_mm_store_ss(&a, _mm_mul_ss(detu, inv_det));
						_mm_store_ss(&b, _mm_mul_ss(detv, inv_det));
						_mm_store_ps(&result.px, _mm_mul_ps(detp, _mm_shuffle_ps(inv_det, inv_det, 0)));
						result.geoN = Vector3(m_preCalcTris[index].nx, m_preCalcTris[index].ny, m_preCalcTris[index].nz).normalized();
					}
					hit = true;
				}
			}
		}
	}
#else
	Vector3 edge[2], tvec, pvec, qvec,rd,ro;				// Implementation of Moller-Trumbore
	rd = Vector3(r.dx,r.dy,r.dz);
	ro = Vector3(r.ox,r.oy,r.oz);
	float det, inv_det;
	ti3 = m_vertexIndices[index];
	edge[0] = m_vertices[ti3.y] - m_vertices[ti3.x];
	edge[1] = m_vertices[ti3.z] - m_vertices[ti3.x];

	pvec = cross(rd, edge[1]);

	det = dot(edge[0], pvec);

#ifdef TEST_CULL
	tvec = ro - m_vertices[ti3.x];
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

	tvec = ro - m_vertices[ti3.x];
	a = dot(tvec, pvec) * inv_det;
	if (a < 0.0 || a > 1.0) return false;

	qvec = cross(tvec, edge[0]);
	b = dot(rd, qvec) * inv_det;
	if (b < 0.0 || a + b > 1.0) return false;

	float newT = dot(edge[1], qvec)*inv_det;
#endif
	if (newT >= tMin && newT < result.t)
	{		
		result.t = newT;
		if (!shadow)
		{
			result.px = ro.x + newT*rd.x; result.py = ro.y + newT*rd.y; result.pz = ro.z + newT*rd.z;
			result.geoN = cross(edge[0], edge[1]).normalized();
		}
		hit = true;
	}
#endif
	if (hit && !shadow)
	{
		ti3 = m_normalIndices[index];
		float c = 1.0f-a-b;
		result.N = Vector3((m_normals[ti3.x]*c+m_normals[ti3.y]*a+m_normals[ti3.z]*b).normalized());
		if (m_texCoordIndices)
		{
			ti3 = m_texCoordIndices[index];
			result.a = m_texCoords[ti3.x].x*c+m_texCoords[ti3.y].x*a+m_texCoords[ti3.z].x*b;
			result.b = m_texCoords[ti3.x].y*c+m_texCoords[ti3.y].y*a+m_texCoords[ti3.z].y*b;
		}
		else
		{
			result.a = a;
			result.b = b;
		}
	}	
	return hit;
}