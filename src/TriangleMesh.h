#ifndef CSE168_TRIANGLE_MESH_H_INCLUDED
#define CSE168_TRIANGLE_MESH_H_INCLUDED

#include "Matrix4x4.h"
#include "Ray.h"

class AABB; 

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

    struct VectorR2
    {
        float x, y;
    };
	
	bool intersect(HitInfo& result, const Ray& r, float tMin = epsilon, float tMax = MIRO_TMAX, u_int index = 0);

	Vector3* m_normals;
	Vector3* m_vertices;
	VectorR2* m_texCoords;

	TupleI3* m_normalIndices;
	TupleI3* m_vertexIndices;
	TupleI3* m_texCoordIndices;
	u_int m_numTris;

	AABB* AABB_PreCalc;				// Precalced AABB's for all triangles, to speed up BVH build.
	void cleanBVHMem();				// Clear the AABB's to save memory (don't need them after build).

protected:
    void loadObj(FILE* fp, const Matrix4x4& ctm);

	ALIGN_SSE struct PrecomputedTriangle		// Used for SSE intersection routine which currently doesn't work...
	{
		float nx, ny, nz, nd;
		float ux, uy, uz, ud;
		float vx, vy, vz, vd;
	};

	PrecomputedTriangle* m_preCalcTris;			// Holds precomputed triangle data for intersection routine.
	bool doPreCalc;
};

#endif // CSE168_TRIANGLE_MESH_H_INCLUDED
