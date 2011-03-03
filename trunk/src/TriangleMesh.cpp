#include "TriangleMesh.h"
#include "Scene.h"
#include "SSE.h"

using namespace std;

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
#ifndef USE_TRI_PACKETS
	if (doPreCalc)
	{
		if (m_preCalcTris) _aligned_free(m_preCalcTris);
		m_preCalcTris = (PrecomputedTriangle*)_aligned_malloc(sizeof(PrecomputedTriangle)*m_numTris, 16);

		Vector3 A, B, C, N, N1, N2, AC, AB;
		VectorR2 Edge1uv, Edge2uv;
		float d, d1, d2;
		m_tangents = (Vector3*)_aligned_malloc(sizeof(Vector3)*m_numTris, 16);
		m_biTangents = (Vector3*)_aligned_malloc(sizeof(Vector3)*m_numTris, 16);

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

			if (m_texCoords)
			{
				ti3 = m_texCoordIndices[i];
				Edge1uv.x = m_texCoords[ti3.y].x - m_texCoords[ti3.x].x;
				Edge1uv.y = m_texCoords[ti3.y].y - m_texCoords[ti3.x].y;
				Edge2uv.x = m_texCoords[ti3.z].x - m_texCoords[ti3.x].x;
				Edge2uv.y = m_texCoords[ti3.z].y - m_texCoords[ti3.x].y;

				float cp = Edge1uv.y * Edge2uv.x - Edge1uv.x * Edge2uv.y;

				if (cp != 0.0f)
				{
					ti3 = m_normalIndices[i];
					float mul = 1.f / cp;
					Vector3 tangent  = ((AB * -Edge2uv.x + AC * Edge1uv.y) * mul).normalized();
					Vector3 normal   = m_normals[ti3.x];
					m_tangents[ti3.x]   = (tangent - normal*dot(normal, tangent)).normalized();
					m_biTangents[ti3.x] = cross(m_tangents[ti3.x], normal);

					normal   = m_normals[ti3.y];
					m_tangents[ti3.y]   = (tangent - normal*dot(normal, tangent)).normalized();
					m_biTangents[ti3.y] = cross(m_tangents[ti3.y], normal);

					normal   = m_normals[ti3.z];
					m_tangents[ti3.z]   = (tangent - normal*dot(normal, tangent)).normalized();
					m_biTangents[ti3.z] = cross(m_tangents[ti3.z], normal);
				}
			}
			
			N1 = cross(AC,N) * Nsq;
			d1 = -dot(N1,A);
			N2 = cross(N,AB) * Nsq;
			d2 = -dot(N2,A);
			m_preCalcTris[i].n[0] = N.x;
			m_preCalcTris[i].n[1] = N.y;
			m_preCalcTris[i].n[2] = N.z;
			m_preCalcTris[i].n[3] = d;
			m_preCalcTris[i].u[0] = N1.x;
			m_preCalcTris[i].u[1] = N1.y;
			m_preCalcTris[i].u[2] = N1.z;
			m_preCalcTris[i].u[3] = d1;
			m_preCalcTris[i].v[0] = N2.x;
			m_preCalcTris[i].v[1] = N2.y;
			m_preCalcTris[i].v[2] = N2.z;
			m_preCalcTris[i].v[3] = d2;
		}
		doPreCalc = false;
	}
#else
	if (doPreCalc && m_texCoords)
	{
		Vector3 A, B, C, AC, AB;
		VectorR2 Edge1uv, Edge2uv;
		m_tangents = (Vector3*)_aligned_malloc(sizeof(Vector3)*m_numTris, 16);
		m_biTangents = (Vector3*)_aligned_malloc(sizeof(Vector3)*m_numTris, 16);

		TriangleMesh::TupleI3 ti3;
		for (u_int i = 0; i < m_numTris; i++)
		{
			ti3 = m_vertexIndices[i];// [m_index];
			A = m_vertices[ti3.x]; //vertex a of triangle
			B = m_vertices[ti3.y]; //vertex b of triangle
			C = m_vertices[ti3.z]; //vertex c of triangle
			AC = C-A;
			AB = B-A;

			ti3 = m_texCoordIndices[i];
			Edge1uv.x = m_texCoords[ti3.y].x - m_texCoords[ti3.x].x;
			Edge1uv.y = m_texCoords[ti3.y].y - m_texCoords[ti3.x].y;
			Edge2uv.x = m_texCoords[ti3.z].x - m_texCoords[ti3.x].x;
			Edge2uv.y = m_texCoords[ti3.z].y - m_texCoords[ti3.x].y;

			float cp = Edge1uv.y * Edge2uv.x - Edge1uv.x * Edge2uv.y;

			if (cp != 0.0f)
			{
				ti3 = m_normalIndices[i];
				float mul = 1.f / cp;
				Vector3 tangent  = ((AB * -Edge2uv.x + AC * Edge1uv.y) * mul).normalized();
				Vector3 normal   = m_normals[ti3.x];
				m_tangents[ti3.x]   = (tangent - normal*dot(normal, tangent)).normalized();
				m_biTangents[ti3.x] = cross(m_tangents[ti3.x], normal);

				normal   = m_normals[ti3.y];
				m_tangents[ti3.y]   = (tangent - normal*dot(normal, tangent)).normalized();
				m_biTangents[ti3.y] = cross(m_tangents[ti3.y], normal);

				normal   = m_normals[ti3.z];
				m_tangents[ti3.z]   = (tangent - normal*dot(normal, tangent)).normalized();
				m_biTangents[ti3.z] = cross(m_tangents[ti3.z], normal);
			}
		}
	}	
	doPreCalc = false;
#endif
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

	bbMin.x = min(A.x, min(B.x, C.x));
	bbMin.y = min(A.y, min(B.y, C.y));
	bbMin.z = min(A.z, min(B.z, C.z));
	bbMax.x = max(A.x, max(B.x, C.x));
	bbMax.y = max(A.y, max(B.y, C.y));
	bbMax.z = max(A.z, max(B.z, C.z));

	*outBox = AABB(bbMin, bbMax);
}

AABB TriangleMesh::getAABB(u_int index)
{
	Vector3 A, B, C, bbMin, bbMax;
	TriangleMesh::TupleI3 ti3;

	ti3 = m_vertexIndices[index];// [m_index];
	A = m_vertices[ti3.x]; //vertex a of triangle
	B = m_vertices[ti3.y]; //vertex b of triangle
	C = m_vertices[ti3.z]; //vertex c of triangle

	bbMin.x = min(A.x, min(B.x, C.x));
	bbMin.y = min(A.y, min(B.y, C.y));
	bbMin.z = min(A.z, min(B.z, C.z));
	bbMax.x = max(A.x, max(B.x, C.x));
	bbMax.y = max(A.y, max(B.y, C.y));
	bbMax.z = max(A.z, max(B.z, C.z));

	AABB outBox = AABB(bbMin, bbMax);
	return outBox;
}