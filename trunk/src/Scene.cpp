#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include <time.h>
#include <omp.h>
#include <Windows.h>

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

int Scene::randsIdx[32] = {0};		// Current index into random number array
float Scene::rands[2097152];	// Array of random numbers
MTRand_int32 Scene::drand(clock());

void Scene::genRands(int threadID)			// Run the random number generator. This way we don't run into
{								// threading problems...
	#pragma omp critical
	for (int i = 65536*threadID; i < 65536*threadID+65536; i++)
	{
		rands[i] = ((float)drand()+0.5) * IntRecip;
	}
}

float Scene::getRand(int threadID)
{
	if (randsIdx[threadID] >= 65536)
	{
		randsIdx[threadID] = 0;
		genRands(threadID);
	}
	return rands[threadID*65536+(randsIdx[threadID]++)];
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

unsigned long int Ray::counter[2048] = {0};
unsigned long int Ray::rayTriangleIntersections[2048] = {0};
unsigned long int BVH::rayBoxIntersections[2048] = {0};

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	clock_t start = clock();		// Added a clock here to time the rendering.

	int remain_horiz		= img->width() % bucket_size;
	int num_horiz_buckets	= img->width() / bucket_size + (remain_horiz > 0);
	int remain_vert			= img->height() % bucket_size;
	int num_vert_buckets	= img->height() / bucket_size + (remain_vert > 0);

	int totalBuckets = num_horiz_buckets*num_vert_buckets;
	bool* doneBuckets = new bool[totalBuckets];
	int shownBuckets = 0;
	bool bucketDone = false;
	for (int i = 0; i < totalBuckets; i++)
	{
		doneBuckets[i] = false;
	}

	static int num_threads = 0;
	HANDLE handle = GetStdHandle(STD_INPUT_HANDLE);
	DWORD events;
	INPUT_RECORD buffer;
	bool abort = false;
	FlushConsoleInputBuffer(handle);

#pragma omp parallel shared (doneBuckets, bucketDone)
	{		
		unsigned int threadID = omp_get_thread_num();
		Ray::counter[threadID] = 0;
		Ray::rayTriangleIntersections[threadID] = 0;
		BVH::rayBoxIntersections[128*threadID] = 0;
		Ray ray(threadID);
		HitInfo hitInfo;

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

							// This is to trap Esc from the console while rendering, to enable aborting the rendering process.
							PeekConsoleInput(handle, &buffer, 1, &events);
							if (events > 0)
							{
								ReadConsoleInput(handle, &buffer, 1, &events);
								if (buffer.Event.KeyEvent.wVirtualKeyCode == VK_ESCAPE)
								{
									abort = true;
									FlushConsoleInputBuffer(handle);
									shownBuckets = totalBuckets;
									#pragma omp flush (abort)
								}
							}
						}
					}
					bucketDone = false;
				} while (bucketDone == false && shownBuckets < totalBuckets);
			}
		}

		#pragma omp for schedule(dynamic) nowait
		for (int bucket = 0; bucket < totalBuckets; bucket++)
		{
			#pragma omp flush (abort)
			if (!abort)
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
					shownBuckets++;
					printf("Rendering Progress: %.3f%%\r", (float)shownBuckets/float(totalBuckets)*100.0f);
					fflush(stdout);
					PeekConsoleInput(handle, &buffer, 1, &events);
					if (events > 0)
					{
						ReadConsoleInput(handle, &buffer, 1, &events);
						if (buffer.Event.KeyEvent.wVirtualKeyCode == VK_ESCAPE)
						{
							abort = true;
							FlushConsoleInputBuffer(handle);
						}
					}
				}
			}
		}
	}
	int counter = 0, rayBoxIntersections = 0, rayTriangleIntersections = 0;
	for (int i = 1; i <= num_threads; i++)
	{
		counter                  += Ray::counter[i*128];
		rayBoxIntersections      += BVH::rayBoxIntersections[i*128];
		rayTriangleIntersections += Ray::rayTriangleIntersections[i*128];
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
	_mm_prefetch((char *)&ray, _MM_HINT_T0);
	_mm_prefetch((char *)&ray.dz4, _MM_HINT_T0);
	Vector3 result = 0;

	hitInfo.t = MIRO_TMAX;
	if (trace(threadID, hitInfo, ray, epsilon))
	{	
		for (int i = 0; i < m_numPaths; i++)
		{
			result += hitInfo.obj->m_material->shade(threadID, ray, hitInfo, *this);
		}
		return result * (1.0f / m_numPaths);
	}
	else
	{
		if (m_envMap != NULL)	// environment map lookup
		{
			result = m_envMap->getLookupXYZ3(ray.d[0], ray.d[1], ray.d[2]) * m_envExposure; 
		}
		else result = m_BGColor;
	}
	return result;
}

int getSum(const int n)
{
	return (int)(n*(n+1)*(2*n+1)*0.16666667f);
}

Vector3 Scene::adaptiveSampleScene(const unsigned int threadID, Camera *cam, Image *img, Ray &ray, HitInfo &hitInfo, int x, int y)
{
	ray = cam->eyeRayAdaptive(threadID, x, y, 0.5f, 0.5f, 0.5f, 0.5f, img->width(), img->height());
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
		Vector3 newResultGamma = Vector3(Image::linear_to_gammaF[int(((newResult.x > 1.f) ? 1.f : newResult.x) * 32767.f)],
										 Image::linear_to_gammaF[int(((newResult.y > 1.f) ? 1.f : newResult.y) * 32767.f)],
										 Image::linear_to_gammaF[int(((newResult.z > 1.f) ? 1.f : newResult.z) * 32767.f)]);
		Vector3 oldResultGamma = Vector3(Image::linear_to_gammaF[int(((shadeResult.x > 1.f) ? 1.f : shadeResult.x) * 32767.f)],
										 Image::linear_to_gammaF[int(((shadeResult.y > 1.f) ? 1.f : shadeResult.y) * 32767.f)],
										 Image::linear_to_gammaF[int(((shadeResult.z > 1.f) ? 1.f : shadeResult.z) * 32767.f)]);

		Vector3 test = oldResultGamma - newResultGamma;
		cutOff = max(fabsf(test.x), max(fabsf(test.y), fabsf(test.z))) < m_noiseThreshold;

		shadeResult = newResult;
		curLevel++;
	}

	return shadeResult;
}

bool Scene::trace(const unsigned int threadID, HitInfo &hitInfo, const Ray& ray, float tMin) const
{
	return m_bvh.intersect(threadID, hitInfo, ray, tMin);
}