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
  m_numPaths = 1;
  m_maxBounces = 10;
  m_noiseThreshold = 0.01f;
  m_minSubdivs = 1;
  m_maxSubdivs = 1;
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

float Scene::getRand()
{
	if (randsIdx >= 990000)
	{
		randsIdx = 0;
		genRands();
	}
	return rands[randsIdx++];
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
        Light* pLight = *lit;
        pLight->preCalc();
    }

    m_bvh.build(&m_objects);
}

unsigned long int Ray::counter[256] = {0};
unsigned long int Ray::rayTriangleIntersections[256] = {0};
unsigned long int BVH::rayBoxIntersections[256] = {0};

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	clock_t start = clock();		// Added a clock here to time the rendering.

	int remain_horiz		= img->width() % bucket_size;
	int num_horiz_buckets	= img->width() / bucket_size + (remain_horiz > 0);
	int remain_vert			= img->height() % bucket_size;
	int num_vert_buckets	= img->height() / bucket_size + (remain_vert > 0);

	static HitInfo hitInfo;
	static Ray ray(1);
	int totalBuckets = num_horiz_buckets*num_vert_buckets;
	bool* doneBuckets = new bool[totalBuckets];
	int shownBuckets = 0;
	bool bucketDone = false;
	for (int i = 0; i < totalBuckets; i++)
	{
		doneBuckets[i] = false;
	}

	static int num_threads = 0;

#pragma omp parallel private(hitInfo, ray) shared (doneBuckets, bucketDone)
	{
		unsigned int threadID = omp_get_thread_num();
		Ray::counter[threadID] = 0;
		Ray::rayTriangleIntersections[threadID] = 0;
		BVH::rayBoxIntersections[threadID] = 0;

		num_threads = omp_get_num_threads();
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

#pragma omp for schedule(dynamic) nowait
		for (int bucket = 0; bucket < totalBuckets; bucket++)
		{
			int bucketX = bucket % num_horiz_buckets;
			int bucketY = bucket / num_horiz_buckets;
			for (int j = bucketY*bucket_size; j < min((bucketY+1)*bucket_size, img->height()); ++j)
			{
				for (int i = bucketX*bucket_size; i < min((bucketX+1)*bucket_size, img->width()); ++i)
				{
					img->setPixel(i, j, adaptiveSampleScene(threadID, cam, img, ray, hitInfo, i, j));
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
	int counter = 0, rayBoxIntersections = 0, rayTriangleIntersections = 0;
	for (int i = 1; i <= num_threads; i++)
	{
		counter                  += Ray::counter[i];
		rayBoxIntersections      += BVH::rayBoxIntersections[i];
		rayTriangleIntersections += Ray::rayTriangleIntersections[i];
	}

	clock_t end = clock();
	printf("Rendering Progress: 100.000%\n");
	debug("done Raytracing!\n");
	printf("Rays cast: %u...\n", counter);
	printf("Ray/AABB intersections: %u...\n", rayBoxIntersections);
	printf("Ray/Triangle intersections: %u...\n", rayTriangleIntersections);
	printf("Rendering time: %.4fs...\n", (end-start)/1000.f);
}

Vector3 Scene::sampleScene(const unsigned int threadID, Ray &ray, HitInfo &hitInfo)
{
	Vector3 result = 0;

	for (int i = 0; i < m_numPaths; i++)
	{
		hitInfo.t = MIRO_TMAX;

		if (trace(threadID, hitInfo, ray, epsilon))
		{					
			result += hitInfo.obj->m_material->shade(threadID, ray, hitInfo, *this);
		}
		else
		{
			if (m_envMap != NULL)	// environment map lookup
			{
				result += m_envMap->getLookupXYZ3(ray.d[0], ray.d[1], ray.d[2]) * m_envExposure; 
			}
			else result += m_BGColor;
		}
	}
	return result * (1.0f / m_numPaths);
}

int getSum(const int n)
{
	return (int)(n*(n+1)*(2*n+1)*0.16666667f);
}

Vector3 Scene::adaptiveSampleScene(const unsigned int threadID, Camera *cam, Image *img, Ray &ray, HitInfo &hitInfo, int x, int y)
{
	ray = cam->eyeRayAdaptive(threadID, x, y, 0.0f, 1.0f, 0.0f, 1.0f, img->width(), img->height());
	Vector3 shadeResult = sampleScene(threadID, ray, hitInfo);

	int curLevel = 2;
	bool cutOff = false;
	while ((curLevel <= m_maxSubdivs && !cutOff) || curLevel <= m_minSubdivs)
	{		
		Vector3 curResult = 0;

		for (int i = 0; i < curLevel; i++)
		{
			for (int j = 0; j < curLevel; j++)
			{
				float offset = 1.0f / (float)curLevel;
				ray = cam->eyeRayAdaptive(threadID, x, y, i*offset, (i+1)*offset, j*offset, (j+1)*offset, img->width(), img->height());
				curResult += sampleScene(threadID, ray, hitInfo);
			}
		}
		float numSamplesPre = getSum(curLevel-1);
		float numSamplesNow = curLevel*curLevel;

		Vector3 newResult = (shadeResult*numSamplesPre + curResult) * (1.0f / (numSamplesPre + numSamplesNow));

		Vector3 test = shadeResult - newResult;
		cutOff = max(fabsf(test.x), max(fabsf(test.y), fabsf(test.z))) < m_noiseThreshold;

		shadeResult = newResult;
		curLevel++;
	}

	return shadeResult;
}

bool Scene::trace(const unsigned int threadID, HitInfo& hitInfo, const Ray& ray, float tMin) const
{
	return m_bvh.intersect(threadID, hitInfo, ray, tMin);
}