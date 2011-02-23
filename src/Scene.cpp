#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include <time.h>
#include <omp.h>

using namespace std;

Scene * g_scene = 0;

Scene::Scene() {
  m_envMap = NULL;
  m_envExposure = 1.0f;
  m_pathTrace = false;
  m_numRays = 1;
  m_maxBounces = 10;
  genRands();
}

int Scene::randsIdx = 0;		// Current index into random number array
float Scene::rands[1000000];	// Array of random numbers
MTRand_int32 Scene::drand(clock());


void Scene::genRands()			// Run the random number generator. This way we don't run into
{								// threading problems...
	#pragma omp critical
	for (int i = 0; i < 1000000; i++)
	{
		rands[i] = ((float)drand()+0.5) * IntRecip;
	}
}

void Scene::openGL(Camera *cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    cam->drawGL();

    // draw objects
    for (size_t i = 0; i < m_objects.size(); ++i)
        m_objects[i]->renderGL();

    glutSwapBuffers();
}

void
Scene::preCalc()
{
    Objects::iterator it;
    for (it = m_objects.begin(); it != m_objects.end(); it++)
    {
        Object* pObject = *it;
        pObject->preCalc();
    }
    Lights::iterator lit;
    for (lit = m_lights.begin(); lit != m_lights.end(); lit++)
    {
        PointLight* pLight = *lit;
        pLight->preCalc();
    }

    m_bvh.build(&m_objects);
}

unsigned long int Ray::counter = 0;
unsigned long int Ray::rayTriangleIntersections = 0;
unsigned long int BVH::rayBoxIntersections = 0;

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	Ray::counter = 0;
	Ray::rayTriangleIntersections = 0;
	BVH::rayBoxIntersections = 0;

	clock_t start = clock();		// Added a clock here to time the rendering.

	int remain_horiz		= img->width() % bucket_size;
	int num_horiz_buckets	= img->width() / bucket_size + (remain_horiz > 0);
	int remain_vert			= img->height() % bucket_size;
	int num_vert_buckets	= img->height() / bucket_size + (remain_vert > 0);

	static HitInfo hitInfo;
	static Ray ray;
	int totalBuckets = num_horiz_buckets*num_vert_buckets;
	bool* doneBuckets = new bool[totalBuckets];
	int shownBuckets = 0;
	bool bucketDone = false;
	for (int i = 0; i < totalBuckets; i++)
	{
		doneBuckets[i] = false;
	}

#pragma omp parallel private(hitInfo, ray) shared (doneBuckets, bucketDone)
	{
		const static int num_threads = omp_get_num_threads();
		if (num_threads > 1)
		{
#pragma omp master
			{
				do		// Draw the buckets as they become ready...
				{
					for (int k = 0; k < totalBuckets; k++)
					{
						if (doneBuckets[k])
						{
							int bucketX = k % num_horiz_buckets;
							int bucketY = k / num_horiz_buckets;
							doneBuckets[k] = false;
							shownBuckets++;
							img->drawBucket(bucketY*bucket_size, min((bucketY+1)*bucket_size, img->height()), bucketX*bucket_size, min((bucketX+1)*bucket_size, img->width()));
							glFinish();
							printf("Rendering Progress: %.3f%%\r", (float)shownBuckets/float(totalBuckets)*100.0f);
							fflush(stdout);
						}
					}
					bucketDone = false;
				} while (bucketDone == false && shownBuckets < totalBuckets);
			}
		}
		float m_numRaysRecip = 1.0f / m_numRays;

#pragma omp for schedule(dynamic) nowait
		for (int bucket = 0; bucket < totalBuckets; bucket++)
		{
			int bucketX = bucket % num_horiz_buckets;
			int bucketY = bucket / num_horiz_buckets;
			for (int j = bucketY*bucket_size; j < min((bucketY+1)*bucket_size, img->height()); ++j)
			{
				for (int i = bucketX*bucket_size; i < min((bucketX+1)*bucket_size, img->width()); ++i)
				{
					hitInfo.t = MIRO_TMAX;
					hitInfo.a = 0.0;
					hitInfo.b = 0.0;
					if (m_pathTrace) {//path tracing
						Vector3 shadeResult = Vector3(0.0f);
						for (int p = 0; p < m_numRays; p++) {
							hitInfo.t = MIRO_TMAX;
							ray = cam->eyeRayRandom(i, j, img->width(), img->height());
							if (trace(hitInfo, ray, epsilon))
							{
								shadeResult += hitInfo.obj->m_material->shade(ray, hitInfo, *this);
							}
							else
							{
								if (m_envMap != NULL) 
								{
									//environment map lookup
									shadeResult += m_envMap->getLookupXYZ3(ray.d[0], ray.d[1], ray.d[2]) * m_envExposure;
								}
								shadeResult += g_scene->getBGColor()*m_envExposure;
							}
						}
						shadeResult *= m_numRaysRecip;
						img->setPixel(i, j, shadeResult);
						
					} else {//no path tracing
						ray = cam->eyeRay(i, j, img->width(), img->height());				
						if (trace(hitInfo, ray, epsilon))
						{					
							img->setPixel(i, j, hitInfo.obj->m_material->shade(ray, hitInfo, *this));
						}
						else
						{
							if (m_envMap != NULL) 
							{
								//environment map lookup
								Vector3 envColor = m_envMap->getLookupXYZ3(ray.d[0], ray.d[1], ray.d[2]);
								envColor *= m_envExposure;
								img->setPixel(i,j,envColor);
							}
							else img->setPixel(i, j, g_scene->getBGColor());
						}
					} //end
				}
			}
			if (num_threads > 1)
			{
				doneBuckets[bucket] = true;
				bucketDone = true;
			}
			else
			{
				img->drawBucket(bucketY*bucket_size, min((bucketY+1)*bucket_size, img->height()), bucketX*bucket_size, min((bucketX+1)*bucket_size, img->width()));
				glFinish();	
			}
		}
	}

	clock_t end = clock();
	printf("Rendering Progress: 100.000%\n");
	debug("done Raytracing!\n");
	printf("Rays cast: %u...\n", Ray::counter);
	printf("Ray/AABB intersections: %u...\n", BVH::rayBoxIntersections);
	printf("Ray/Triangle intersections: %u...\n", Ray::rayTriangleIntersections);
	printf("Rendering time: %.4fs...\n", (end-start)/1000.f);
}

bool Scene::trace(HitInfo& minHit, const Ray& ray, float tMin) const
{
	return m_bvh.intersect(minHit, ray, tMin);
}