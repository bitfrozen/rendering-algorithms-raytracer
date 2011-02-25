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

	const Vector3 getPointOnLight() const { return m_position; }
	
	void preCalc() {} // use this if you need to

protected:
    Vector3 m_position;
};

#endif // CSE168_POINTLIGHT_H_INCLUDED
