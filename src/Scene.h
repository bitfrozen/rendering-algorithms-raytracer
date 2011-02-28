#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Light.h"
#include "Miro.h"
#include "Object.h"
#include "BVH.h"
#include "MTRand.h"

class Camera;
class Image;

class Scene
{
public:
	Scene();
    void addObject(Object* pObj)        {m_objects.push_back(pObj);}
    const Objects* objects() const      {return &m_objects;}

    void addLight(Light* pObj)			{m_lights.push_back(pObj);}
    const Lights* lights() const        {return &m_lights;}

	void setEnvMap(Texture* envMap)		{m_envMap = envMap; }
	void setEnvExposure(float exposure)	{m_envExposure = exposure;}
	const float getEnvExposure() const  {return m_envExposure;}
	const Texture* getEnvMap() const	{return m_envMap;}

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    bool trace(const unsigned int threadID, HitInfo& hitInfo, const Ray& ray, float tMin = epsilon) const;

	Vector3 adaptiveSampleScene(const unsigned int threadID, Camera *cam, Image *img, Ray &ray, HitInfo &hitInfo, int i, int j);
	Vector3 sampleScene(const unsigned int threadID, Ray &ray, HitInfo &hitInfo);

	void setBGColor(Vector3& color)		{m_BGColor = color;}
	const Vector3& getBGColor()	const   {return m_BGColor;}

	void setPathTrace(bool pt)		{ m_pathTrace = pt;}

	void setMinSubdivs(int r)		{ m_minSubdivs = r;}
	int minSubdivs()				{ return m_minSubdivs;}

	void setMaxSubdivs(int r)		{ m_maxSubdivs = r;}
	int maxSubdivs()				{ return m_maxSubdivs;}

	void setMaxBounces(int mb)		{ m_maxBounces = mb;}
	int maxBounces()				{ return m_maxBounces;}

	void setNumPaths(int p)			{ m_numPaths = p;}
	int numPaths()					{ return m_numPaths;}

	void setNoise(const float n)	{ m_noiseThreshold = n;}
	const float noise() const		{ return m_noiseThreshold;}
	
	bool m_pathTrace;
	int m_numPaths;
	int m_minSubdivs;
	int m_maxSubdivs;
	int m_maxBounces;

	static int randsIdx;
	static float rands[1000000];
	static MTRand_int32 drand;
	static void genRands();
	static float getRand();

protected:
    Objects m_objects;
    BVH m_bvh;
    Lights m_lights;
	Vector3 m_BGColor;
	Texture* m_envMap;
	float m_envExposure;
	float m_noiseThreshold;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
