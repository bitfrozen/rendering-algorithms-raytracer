#ifndef CSE168_RECT_LIGHT_H_INCLUDED
#define CSE168_RECT_LIGHT_H_INCLUDED

#include "Light.h"
#include "Scene.h"

// Implements a direct sampled arbitrary parallelogram light
class RectangleLight : public Light
{
public:
	RectangleLight() : Light() {m_v1 = m_v2 = m_v3 = Vector3(0); setType(RECTANGLE_LIGHT);}

    void setVertices(const Vector3& v1, const Vector3& v2, const Vector3& v3);
	void setPower(float f);
    
	const Vector3 sampleLight(const unsigned int threadID, const Vector3 &from, const Vector3 &normal, const Scene &scene, const Vector3 &rVec, float &outSpec) const;

protected:
    Vector3 m_v1, m_v2, m_v3;
};

#endif