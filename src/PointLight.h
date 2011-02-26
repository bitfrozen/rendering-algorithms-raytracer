#ifndef CSE168_POINTLIGHT_H_INCLUDED
#define CSE168_POINTLIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Light.h"

class PointLight : public Light
{
public:
    void setPosition(const Vector3& v)  {m_position = v;}
    
    const Vector3& position() const     {return m_position;}

	const Vector3 sampleLight(const Vector3 &from, const Vector3 &normal, const Scene &scene, const Vector3 &rVec, float &outSpec) const;

protected:
    Vector3 m_position;
};

#endif // CSE168_POINTLIGHT_H_INCLUDED