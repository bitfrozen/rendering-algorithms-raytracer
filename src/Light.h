#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"

class Light
{
public:
    void setColor(const Vector3& v)     {m_color = v;}
    void setWattage(float f)            {m_wattage = f;}
    
    float wattage() const               {return m_wattage;}
    const Vector3& color() const        {return m_color;}

	const Vector3 getPointOnLight() const {}

    virtual void preCalc() {}

protected:
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED