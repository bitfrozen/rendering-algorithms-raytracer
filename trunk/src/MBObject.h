#ifndef CSE168_MBOBJECT_H_INCLUDED
#define CSE168_MBOBJECT_H_INCLUDED

#include "Object.h"
#include "Miro.h"
#include "Material.h"
#include <vector>

using namespace std;

ALIGN_SSE class MBObject : public Object
{
public:
	MBObject(Material* mat = NULL, TriangleMesh* m = NULL, TriangleMesh* m2 = NULL, u_int i = 0);
	~MBObject();

	void setMeshT2(TriangleMesh* m) {m_mesh_t2 = m;}

	void preCalc();

	const bool intersect(const unsigned int threadID, HitInfo& result, const Ray& ray, const float tMin = epsilon);

	void getAABB(AABB* outBox);
	AABB getAABB();

	TriangleMesh* m_mesh_t2;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_MBOBJECT_H_INCLUDED