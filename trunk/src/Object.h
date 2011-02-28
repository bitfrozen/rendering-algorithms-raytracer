#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include "Miro.h"
#include "Material.h"
#include <vector>

using namespace std;

ALIGN_SSE class AABB
{
public:
	inline AABB() : bbMin(MIRO_TMAX), bbMax(-MIRO_TMAX) {}
	inline AABB(const Vector3& bbMin, const Vector3& bbMax) : bbMin(bbMin), bbMax(bbMax) {
	}
	inline AABB(const AABB& bb1, const AABB& bb2) {
		bbMin.x = min(bb1.bbMin.x, bb2.bbMin.x);
		bbMin.y = min(bb1.bbMin.y, bb2.bbMin.y);
		bbMin.z = min(bb1.bbMin.z, bb2.bbMin.z);
		bbMax.x = max(bb1.bbMax.x, bb2.bbMax.x);
		bbMax.y = max(bb1.bbMax.y, bb2.bbMax.y);
		bbMax.z = max(bb1.bbMax.z, bb2.bbMax.z);
	}
	inline void grow(const Vector3& newPt) {
		bbMin.x = min(bbMin.x, newPt.x);
		bbMin.y = min(bbMin.y, newPt.y);
		bbMin.z = min(bbMin.z, newPt.z);
		bbMax.x = max(bbMax.x, newPt.x);
		bbMax.y = max(bbMax.y, newPt.y);
		bbMax.z = max(bbMax.z, newPt.z);
	}
	inline float getArea() {
		return 2*((bbMax.x-bbMin.x + bbMax.z-bbMin.z)*(bbMax.y-bbMin.y) + (bbMax.x-bbMin.x)*(bbMax.z-bbMin.z));
	}
	inline Vector3 getCentroid() {
		return 0.5f * (bbMin+bbMax);
	}

#ifndef NO_SSE
	ALIGN_SSE Vector3 bbMin;
	ALIGN_SSE Vector3 bbMax;
#else
	Vector3 bbMin, bbMax;
#endif
};

enum objectType_t {OBJECT, MB_OBJECT};

ALIGN_SSE class Object
{
public:
    Object(Material* mat = NULL, TriangleMesh* m = NULL, u_int i = 0);
    ~Object();

    void setMaterial(const Material* m) {m_material = m;}
	void setIndex(u_int i) {m_index = i;}
	void setMesh(TriangleMesh* m) {m_mesh = m;}

    void renderGL();
    virtual void preCalc();

    const virtual bool intersect(const unsigned int threadID, HitInfo& result, const Ray& ray, const float tMin = epsilon);

	static int sortByXComponent(const void* s1, const void* s2);	// Sorting functions for use with qsort
	static int sortByYComponent(const void* s1, const void* s2);
	static int sortByZComponent(const void* s1, const void* s2);
	static int sortByArea(const void* s1, const void* s2);

	virtual void getAABB(AABB* outBox);
	virtual AABB getAABB();

	const Material* m_material;
	TriangleMesh* m_mesh;
	u_int m_index;
	objectType_t m_objectType;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED