#ifndef CSE168_TRIANGLE_H_INCLUDED
#define CSE168_TRIANGLE_H_INCLUDED

#include "Object.h"

/*
    The Triangle class stores a pointer to a mesh and an index into its
    triangle array. The mesh stores all data needed by this Triangle.
*/
class Triangle : public Object
{
public:
    Triangle(TriangleMesh* m = 0, u_int i = 0);
    virtual ~Triangle();

	virtual void preCalc();
    void setIndex(u_int i) {m_index = i;}
    void setMesh(TriangleMesh* m) {m_mesh = m;}

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = epsilon, float tMax = MIRO_TMAX);

	virtual void getAABB(AABB& bBox);
	virtual void getCentroid(Vector3& centroid);
	void cleanBVHMem();

protected:
    TriangleMesh* m_mesh;
    u_int m_index;
};

#endif // CSE168_TRIANGLE_H_INCLUDED
