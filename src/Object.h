#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
#include "Material.h"

ALIGN_SSE class AABB
{
public:
	AABB() : bbMin(1), bbMax(-1), area(MIRO_TMAX), centroid(0) {}
	AABB(const Vector3& bbMin, const Vector3& bbMax) : bbMin(bbMin), bbMax(bbMax) {
		area = 2*((bbMax.x-bbMin.x + bbMax.z-bbMin.z)*(bbMax.y-bbMin.y) + (bbMax.x-bbMin.x)*(bbMax.z-bbMin.z));
		centroid =  0.5f * (bbMin+bbMax);
	}
	AABB(const AABB& bb1, const AABB& bb2) {
		bbMin.x = std::min(bb1.bbMin.x, bb2.bbMin.x);
		bbMin.y = std::min(bb1.bbMin.y, bb2.bbMin.y);
		bbMin.z = std::min(bb1.bbMin.z, bb2.bbMin.z);
		bbMax.x = std::max(bb1.bbMax.x, bb2.bbMax.x);
		bbMax.y = std::max(bb1.bbMax.y, bb2.bbMax.y);
		bbMax.z = std::max(bb1.bbMax.z, bb2.bbMax.z);
		area = 2*((bbMax.x-bbMin.x + bbMax.z-bbMin.z)*(bbMax.y-bbMin.y) + (bbMax.x-bbMin.x)*(bbMax.z-bbMin.z));
		centroid =  0.5f * (bbMin+bbMax);
	}

#ifndef NO_SSE
	ALIGN_SSE Vector3 bbMin;
	ALIGN_SSE Vector3 bbMax;
#else
	Vector3 bbMin, bbMax;
#endif
	float area;
	Vector3 centroid;
};

class Object
{
public:
    Object() {m_bBox = new AABB;}
    virtual ~Object() {delete m_bBox;}

    void setMaterial(const Material* m) {m_material = m;}
	const Material* getMaterial() const	{return m_material;}

    virtual void renderGL() {}
    virtual void preCalc() {}

    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX) = 0;

	static int sortByXComponent(const void* s1, const void* s2);	// Sorting functions for use with qsort (defined in BVH.cpp)
	static int sortByYComponent(const void* s1, const void* s2);
	static int sortByZComponent(const void* s1, const void* s2);
	AABB* m_bBox;

protected:
    const Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
