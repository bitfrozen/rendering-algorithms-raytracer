#ifndef CSE168_TRIANGLE_MESH_H_INCLUDED
#define CSE168_TRIANGLE_MESH_H_INCLUDED

#include "Matrix4x4.h"
#include "Ray.h"
#include "Object.h"

class TriangleMesh
{
public:
    TriangleMesh();
    ~TriangleMesh();

    // load from an OBJ file
    bool load(char* file, const Matrix4x4& ctm = Matrix4x4());

    // for single triangles
    void createSingleTriangle();
	void preCalc();
    inline void setV1(const Vector3& v) {m_vertices[0] = v;}
    inline void setV2(const Vector3& v) {m_vertices[1] = v;}
    inline void setV3(const Vector3& v) {m_vertices[2] = v;}
    inline void setN1(const Vector3& n) {m_normals[0] = n;}
    inline void setN2(const Vector3& n) {m_normals[1] = n;}
    inline void setN3(const Vector3& n) {m_normals[2] = n;}

    struct TupleI3
    {
        u_int x, y, z;
    };

    ALIGN_SSE struct VectorR2
    {
        float x, y;
    };

	void getAABB(u_int index, AABB* outBox);
	AABB getAABB(u_int index);

	ALIGN_SSE Vector3* m_normals;
	ALIGN_SSE Vector3* m_vertices;
	ALIGN_SSE VectorR2* m_texCoords;

	TupleI3* m_normalIndices;
	TupleI3* m_vertexIndices;
	TupleI3* m_texCoordIndices;
	u_int m_numTris;
	bool doPreCalc;
	
	ALIGN_SSE struct PrecomputedTriangle		// Used for SSE intersection routine
	{
		union {float n[4]; __m128 _n;};
		union {float u[4]; __m128 _u;};
		union {float v[4]; __m128 _v;};
	};
	PrecomputedTriangle* m_preCalcTris;			// Holds precomputed triangle data for intersection routine.

protected:
    void loadObj(FILE* fp, const Matrix4x4& ctm);
};

#endif // CSE168_TRIANGLE_MESH_H_INCLUDED
