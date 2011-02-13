#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include <vector>

#define GET_NUMCHILD(a) (a)>>1

class BVH_Node
{
public:
	AABB* bBox;
	union {
		u_int numChildren;
		u_int isLeaf;
	};
	union {
		BVH_Node* Children;
		Object** objs;
	};
	bool intersect(HitInfo& result, const Ray& ray, float tMin = epsilon, float tMax = MIRO_TMAX) const;
	void build(Object** objs, AABB* preCalcAABB, u_int numObjs, float* leftArea, float* rightArea, int* binIds);
	void partitionSweep(Object** objs, AABB* preCalcAABB, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, int* binIds);
};

class BVH
{
public:
    void build(Objects* objs);

    bool intersect(HitInfo& result, const Ray& ray,
                   float tMin = epsilon, float tMax = MIRO_TMAX) const;
protected:
    Objects* m_objects;
	BVH_Node m_baseNode;
};

#endif // CSE168_BVH_H_INCLUDED
