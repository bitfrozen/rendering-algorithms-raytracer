#ifndef CSE168_SPHERE_H_INCLUDED
#define CSE168_SPHERE_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Sphere : public Object
{
public:
	Sphere();
    virtual ~Sphere();

	void setCenter(const Vector3& v)    {m_center = v; m_AABB = AABB(m_center-m_radius, m_center+m_radius);}
	void setRadius(const float f)       {m_radius = f; m_AABB = AABB(m_center-m_radius, m_center+m_radius);}

    const Vector3& center() const       {return m_center;}
    float radius() const                {return m_radius;}
	virtual void getAABB(AABB* outBox)	{outBox = &m_AABB;}
	virtual AABB getAABB()				{return m_AABB;}

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

protected:
    Vector3 m_center;
    float m_radius;
	AABB m_AABB;
};

#endif // CSE168_SPHERE_H_INCLUDED
