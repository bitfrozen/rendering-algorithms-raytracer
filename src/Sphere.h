#ifndef CSE168_SPHERE_H_INCLUDED
#define CSE168_SPHERE_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Sphere : public Object
{
public:
	Sphere();
    virtual ~Sphere();

    void setCenter(const Vector3& v)    {m_center = v; *bBox = AABB(m_center-m_radius, m_center+m_radius);}
    void setRadius(const float f)       {m_radius = f; *bBox = AABB(m_center-m_radius, m_center+m_radius);}

    const Vector3& center() const       {return m_center;}
    float radius() const                {return m_radius;}

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

	virtual void getAABB(AABB& box) {box = *bBox;}
	virtual void getCentroid(Vector3& centroid) {centroid = m_center;}

protected:
    Vector3 m_center;
    float m_radius;
	AABB* bBox;
};

#endif // CSE168_SPHERE_H_INCLUDED
