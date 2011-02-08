#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include <time.h>
#include <omp.h>

Scene * g_scene = 0;

Scene::Scene() {
  m_envMap = NULL;
  m_envExposure = 1.0f;
}

void
Scene::openGL(Camera *cam)
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

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	clock_t start = clock();		// Added a clock here to time the rendering.
    Ray ray;
    HitInfo hitInfo;
    Vector3 shadeResult;

#pragma omp parallel private(ray, hitInfo, shadeResult)
	{
		// loop over all pixels in the image
		for (int j = 0; j < img->height(); ++j)
		{
#pragma omp for schedule(dynamic, 4)
			for (int i = 0; i < img->width(); ++i)
			{
				hitInfo.t = MIRO_TMAX;
				ray = cam->eyeRay(i, j, img->width(), img->height());
				if (trace(hitInfo, ray, 0.0001, MIRO_TMAX))
				{
					shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
					img->setPixel(i, j, shadeResult);
				} 
				else 
				{
					if (m_envMap != NULL) 
					{
						//environment map lookup
						Vector3 rayD = Vector3(ray.dx, ray.dy, ray.dz);
						float theta = atan2(rayD.z, rayD.x) + PI;
						float phi = acos(rayD.y);
						float u = theta * 0.5 * piRecip;
						float v = 1.0 - (phi * piRecip);

						Vector3 envColor = m_envMap->getLookup3(u,v);
						envColor /= m_envExposure;
						img->setPixel(i,j,envColor);
					}
					else img->setPixel(i, j, g_scene->getBGColor());
				}
			}
			img->drawScanline(j);
			glFinish();
			printf("Rendering Progress: %.3f%%\r", j/float(img->height())*100.0f);
			fflush(stdout);
		}
	}
	clock_t end = clock();
	printf("Rendering Progress: 100.000%\n");
	debug("done Raytracing!\n");
	printf("Rendering time: %.4fs...\n", (end-start)/1000.f);
}

bool
Scene::trace(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    return m_bvh.intersect(minHit, ray, tMin, tMax);
}
