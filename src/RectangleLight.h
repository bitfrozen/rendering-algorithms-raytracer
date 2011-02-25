#ifndef CSE168_RECT_LIGHT_H_INCLUDED
#define CSE168_RECT_LIGHT_H_INCLUDED

#include "Light.h"
#include "Scene.h"

class RectangleLight : public Light
{
public:
    void setVertices(const Vector3& v1, const Vector3& v2, const Vector3& v3)  {m_v1 = v1; m_v2 = v2; m_v3 = v3; }
    
	const Vector3 getPointOnLight() const;
	
	void preCalc() {} // use this if you need to

protected:
    Vector3 m_v1;
	Vector3 m_v2;
	Vector3 m_v3;

};

#endif