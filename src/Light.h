#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Ray.h"
#include "Scene.h"

enum lightType_t {RECTANGLE_LIGHT, POINT_LIGHT, DOME_LIGHT};

class Light
{
public:
	Light() {
		m_color = Vector3(0); m_power = 0; m_numSamples = 1; 
		m_lightType = POINT_LIGHT; m_castShadows = true; m_fastShadows = true;
		m_noiseThreshold = epsilon;
	}
    void setColor(const Vector3& v)     {m_color = v;}
    virtual void setPower(float f)      {m_power = f;}
	void setSamples(int n)				{m_numSamples = n;}
	void setType(lightType_t newType)	{m_lightType = newType;}
	void setCastShadows(bool c)			{m_castShadows = c;}
	void setFastShadows(bool c)			{m_fastShadows = c;}
	void setNoiseThreshold(float t)		{m_noiseThreshold = t;}
    
    virtual const float power() const   {return m_power;}
    const Vector3& color() const        {return m_color;}
	int samples() const					{return m_numSamples;}
	lightType_t type() const			{return m_lightType;}
	bool castShadows() const			{return m_castShadows;}
	bool fastShadows() const			{return m_fastShadows;}
	float noiseThreshold() const		{return m_noiseThreshold;}

	const virtual Vector3 sampleLight(const unsigned int threadID, const Vector3 &from, const Vector3 &normal, const float time, const Scene &scene, const Vector3 &rVec, float &outSpec) const = 0;

    virtual void preCalc()				{}

protected:
    Vector3 m_color;
    float m_power;
	int m_numSamples;
	lightType_t m_lightType;
	bool m_castShadows;
	bool m_fastShadows;
	float m_noiseThreshold;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED