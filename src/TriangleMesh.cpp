#include "TriangleMesh.h"
#include "Triangle.h"
#include "Scene.h"
#include <xmmintrin.h>
#include <smmintrin.h>

#define NO_SSE

TriangleMesh::TriangleMesh() :
    m_normals(0),
    m_vertices(0),
    m_texCoords(0),
    m_normalIndices(0),
    m_vertexIndices(0),
    m_texCoordIndices(0),
	m_preCalcTris(0)
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
	m_preCalcTris = new PrecomputedTriangle[m_numTris];
	Vector3 A,B,C,N,N1,N2,AC,AB;
	float d, d1, d2;
	TriangleMesh::TupleI3 ti3;
	for (unsigned int i = 0; i < m_numTris; i++)
	{
		ti3 = m_vertexIndices[i];// [m_index];
		A = m_vertices[ti3.x]; //vertex a of triangle
		B = m_vertices[ti3.y]; //vertex b of triangle
		C = m_vertices[ti3.z]; //vertex c of triangle
		AC = C-A;
		AB = B-A;
		N = cross(AB,AC);
		d = dot(A,N);
		float Nsq = 1 / N.length2();
		
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
}

bool TriangleMesh::intersect(HitInfo& result, const Ray& r, float tMin /* = 0.0f */, float tMax /* = MIRO_TMAX */)
{
	Hit2 theHit;
	theHit.t = MIRO_TMAX;
	Ray2 theRay;
	theRay.ox = r.o.x;
	theRay.oy = r.o.y;
	theRay.oz = r.o.z;
	theRay.ow = 1.0;
	theRay.dx = r.d.x;
	theRay.dy = r.d.y;
	theRay.dz = r.d.z;
	theRay.dw = 0.0;
	bool hit = false;
	TriangleMesh::TupleI3 ti3;
	HitInfo temp;
	for (unsigned int i = 0; i < m_numTris; ++i)
	{
#ifndef NO_SSE
		if (singleIntersect(theHit, theRay, tMin, tMax, i))
		{
			hit = true;
			if (theHit.t < result.t)
			{
				ti3 = m_normalIndices[i];// [m_index]; 
				result.t = theHit.t;
				result.P = Vector3(theHit.px, theHit.py, theHit.pz);
				result.N = (m_normals[ti3.x]*(1-theHit.u-theHit.v)+m_normals[ti3.y]*theHit.u+m_normals[ti3.z]*theHit.v).normalized();
			}
		}
#else
		if (singleIntersect(temp, r, tMin, tMax, i))
		{
			hit = true;
			if (temp.t < result.t)
			{
				result = temp;
			}
		}
#endif
	}
	return hit;
}

bool TriangleMesh::singleIntersect(Hit2& result, const Ray2& r,float tMin, float tMax, int index)
{
	__declspec(align(16)) const float arr[4] = {-1,-1,-1,1};
	const __m128 int_coef = _mm_load_ps(arr);

	const __m128 o = _mm_load_ps(&r.ox);
	const __m128 d = _mm_load_ps(&r.dx);
	const __m128 n = _mm_load_ps(&m_preCalcTris[index].nx);

	const __m128 det = _mm_dp_ps(n, d, 0x7f);
	const __m128 dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
	const __m128 oldt = _mm_load_ss(&result.t);

	if ((_mm_movemask_ps(_mm_xor_ps(dett,_mm_sub_ss(_mm_mul_ss(oldt, det), dett)))&1) == 0)
	{
		const __m128 detp = _mm_add_ps(_mm_mul_ps(o, det), _mm_mul_ps(dett, d));
		const __m128 detu = _mm_dp_ps(detp, _mm_load_ps(&m_preCalcTris[index].ux), 0xf1);
		if ((_mm_movemask_ps(_mm_xor_ps(detu, _mm_sub_ss(det, detu)))&1) == 0)
		{
			const __m128 detv = _mm_dp_ps(detp, _mm_load_ps(&m_preCalcTris[index].vx), 0xf1);

			if ((_mm_movemask_ps(_mm_xor_ps(detv, _mm_sub_ss(det, _mm_add_ss(detu, detv))))&1) == 0)
			{
				const __m128 inv_det = _mm_rcp_ss(det);
				_mm_store_ss(&result.t, _mm_mul_ss(dett, inv_det));
				_mm_store_ss(&result.u, _mm_mul_ss(detu, inv_det));
				_mm_store_ss(&result.v, _mm_mul_ss(detv, inv_det));
				_mm_store_ps(&result.px, _mm_mul_ps(detp, _mm_shuffle_ps(inv_det, inv_det, 0)));
				return true;
			}
		}
	}
	return false;
}

bool
TriangleMesh::singleIntersect(HitInfo& result, const Ray& r, float tMin, float tMax, int index)
{
	Vector3 edge1, edge2, tvec, pvec, qvec;
	float det, inv_det, u, v;
	TriangleMesh::TupleI3 ti3;
	ti3 = m_vertexIndices[index];
	edge1 = m_vertices[ti3.y] - m_vertices[ti3.x];
	edge2 = m_vertices[ti3.z] - m_vertices[ti3.x];

	pvec = cross(r.d, edge2);

	det = dot(edge1, pvec);

	tvec = r.o - m_vertices[ti3.x];
	u = dot(tvec, pvec);
	if (u < 0.0 || u > det) return false;

	qvec = cross(tvec, edge1);
	v = dot(r.d, qvec);
	if (v < 0.0 || u + v > det) return false;

	result.t = dot(edge2, qvec);
	inv_det = 1.0 / det;
	result.t *= inv_det;
	u *= inv_det;
	v *= inv_det;

	ti3 = m_normalIndices[index];
	result.N = Vector3((m_normals[ti3.x]*(1-u-v)+m_normals[ti3.y]*u+m_normals[ti3.z]*v).normalized());
	return true;
}