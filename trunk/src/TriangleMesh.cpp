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
	delete [] m_preCalcTris;
}

void TriangleMesh::preCalc()
{
#ifndef NO_SSE
	if (doPreCalc)			// Precalc a few values to speed up BVH build. The rest of the data was for SSE but it's not working, so...
	{
		delete[] m_preCalcTris;
		m_preCalcTris = (PrecomputedTriangle*)_aligned_malloc(sizeof(PrecomputedTriangle)*m_numTris, 16); //new PrecomputedTriangle[m_numTris];

		Vector3 A, B, C, N,N1,N2,AC,AB;
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

bool TriangleMesh::intersect(HitInfo& result, const Ray& r, float tMin, float tMax, u_int index)
{
	bool hit = false;
	TriangleMesh::TupleI3 ti3;
#ifndef NO_SSE													// Implementation of Havel-Herout(2010), works fine for diffuse shading but breaks up for reflection/refraction
	ALIGN_SSE float newT, u, v;
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
				_mm_store_ss(&u, _mm_mul_ss(detu, inv_det));
				_mm_store_ss(&v, _mm_mul_ss(detv, inv_det));
				if (newT >= tMin && newT < result.t && u >= 0.0f && v >= 0.0f && u+v <= 1.0f)
				{
					result.t = newT;
					_mm_store_ps(&result.px, _mm_mul_ps(detp, _mm_shuffle_ps(inv_det, inv_det, 0)));
					hit = true;
				}
			}
		}
	}
	
	if (hit)
	{
		ti3 = m_normalIndices[index];
		result.N = (m_normals[ti3.x]*(1-u-v)+m_normals[ti3.y]*u+m_normals[ti3.z]*v).normalized();
		result.geoN = Vector3(m_preCalcTris[index].nx, m_preCalcTris[index].ny, m_preCalcTris[index].nz).normalized();
		if (m_texCoordIndices)
		{
			ti3 = m_texCoordIndices[index];
			result.a = m_texCoords[ti3.x].x*(1-u-v)+m_texCoords[ti3.y].x*u+m_texCoords[ti3.z].x*v;
			result.b = m_texCoords[ti3.x].y*(1-u-v)+m_texCoords[ti3.y].y*u+m_texCoords[ti3.z].y*v;
		}
		else
		{
			result.a = u;
			result.b = v;
		}
	}
#else
	Vector3 edge[2], tvec, pvec, qvec,rd,ro;				// Implementation of Moller-Trumbore
	rd = Vector3(r.dx,r.dy,r.dz);
	ro = Vector3(r.ox,r.oy,r.oz);
	float det, inv_det, u, v;
	ti3 = m_vertexIndices[index];
	edge[0] = m_vertices[ti3.y] - m_vertices[ti3.x];
	edge[1] = m_vertices[ti3.z] - m_vertices[ti3.x];

	pvec = cross(rd, edge[1]);

	det = dot(edge[0], pvec);

#ifdef TEST_CULL
	tvec = ro - m_vertices[ti3.x];
	u = dot(tvec, pvec);	
	if (u < 0.0 || u > det) return false;

	qvec = cross(tvec, edge[0]);
	v = dot(rd, qvec);
	if (v < 0.0 || u + v > det) return false;

	inv_det = 1.0 / det;
	float newT = dot(edge[1], qvec)*inv_det;
	u *= inv_det;
	v *= inv_det;
#else
	inv_det = 1.0 / det;

	tvec = ro - m_vertices[ti3.x];
	u = dot(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0) return false;

	qvec = cross(tvec, edge[0]);
	v = dot(rd, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0) return false;

	float newT = dot(edge[1], qvec)*inv_det;
#endif
	if (newT >= tMin && newT < result.t)
	{		
		result.t = newT;

		result.geoN = cross(edge[0], edge[1]).normalized();
		ti3 = m_normalIndices[index];
		result.N = Vector3((m_normals[ti3.x]*(1-u-v)+m_normals[ti3.y]*u+m_normals[ti3.z]*v).normalized());
		if (m_texCoordIndices)
		{
			ti3 = m_texCoordIndices[index];
			result.a = m_texCoords[ti3.x].x*(1-u-v)+m_texCoords[ti3.y].x*u+m_texCoords[ti3.z].x*v;
			result.b = m_texCoords[ti3.x].y*(1-u-v)+m_texCoords[ti3.y].y*u+m_texCoords[ti3.z].y*v;
		}
		else
		{
			result.a = u;
			result.b = v;
		}
		Vector3 P = ro + result.t*rd;
		result.px = P.x; result.py = P.y; result.pz = P.z;
		hit = true;
	}
#endif
	return hit;
}