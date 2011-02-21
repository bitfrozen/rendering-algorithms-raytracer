#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "BVH.h"

class Camera;
class Image;

class Scene
{
public:
	Scene();
    void addObject(Object* pObj)        {m_objects.push_back(pObj);}
    const Objects* objects() const      {return &m_objects;}

    void addLight(PointLight* pObj)     {m_lights.push_back(pObj);}
    const Lights* lights() const        {return &m_lights;}

	void setEnvMap(Texture* envMap)		{m_envMap = envMap; }
	void setEnvExposure(float exposure)	{m_envExposure = exposure;}
	float getEnvExposure()				{return m_envExposure;}

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray, float tMin = 0.0f) const;

	void setBGColor(Vector3& color)		{m_BGColor = color;}
	Vector3& getBGColor()				{return m_BGColor;}

	Texture* getEnvMap()				{return m_envMap;}

	void setPathTrace(bool pt) { m_pathTrace = pt; }
	void setNumRays(int r) { m_numRays = r; }
	void setMaxBounces(int mb) { m_maxBounces = mb; }
	
	bool m_pathTrace;
	int m_numRays;
	int m_maxBounces;

protected:
    Objects m_objects;
    BVH m_bvh;
    Lights m_lights;
	Vector3 m_BGColor;
	Texture* m_envMap;
	float m_envExposure;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
