#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include <vector>

class BVH_Node
{
public:
	bool intersect(HitInfo& result, const Ray& ray, float tMin = 0.0f, float tMax = MIRO_TMAX) const;
	void build(Object* obj, int first = 0, int last = -100);
	Vector3 vMin, vMax;
	std::vector<BVH_Node> Children;
	bool isLeaf;
	Object* obj;
};

class BVH
{
public:
    void build(Objects * objs);

    bool intersect(HitInfo& result, const Ray& ray,
                   float tMin = 0.0f, float tMax = MIRO_TMAX) const;
protected:
    Objects * m_objects;
	BVH_Node m_baseNode;
};

#endif // CSE168_BVH_H_INCLUDED
