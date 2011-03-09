#include "Scene.h"
#include "Camera.h"
#include "DomeLight.h"
#include "SSE.h"

using namespace std;

void DomeLight::setTexture(Texture* t)
{
	m_lightMap = t;
	if (cosTableU) delete[] cosTableU;
	if (sinTableU) delete[] sinTableU;
	if (cosTableV) delete[] cosTableV;
	if (sinTableV) delete[] sinTableV;
	if (vDistribs)
	{
		for (int i = 0; i < uDistrib->count; ++i)
		{
			if (vDistribs[i]) delete vDistribs[i];
		}
		delete vDistribs;
	}
	if (uDistrib) delete uDistrib;

	int nu = t->getWidth(), nv = t->getHeight();
	float *img = new float[nu*nv];
	for (int u = 0; u < nu; ++u)
	{
		float up = (float)u / (float)nu;
		for (int v = 0; v < nv; ++v)
		{
			float vp = (float)v / (float)nv;
			Vector3 light = m_lightMap->getLookup3(up, vp);
			img[v+u*nv] = light.average();
			if (img[v+u*nv] > 10000)
			{
				int tmp = 0;
			}
		}
	}

	float *func = new float[max(nu, nv)];
	float *sinVals = new float[nv];
	for (int i = 0; i < nv; ++i)
	{
		sinVals[i] = sin(PI * float(i+.5)/float(nv));
	}

	vDistribs = new Distribution1D*[nu];
	for (int u = 0; u < nu; ++u)
	{
		for (int v = 0; v < nv; ++v)
			func[v] = img[u*nv+v] * sinVals[v];
		vDistribs[u] = new Distribution1D(func, nv);
	}

	for (int u = 0; u < nu; ++u)
		func[u] = vDistribs[u]->funcInt;
	uDistrib = new Distribution1D(func, nu);

	cosTableU = new float[nu+1];
	sinTableU = new float[nu+1];
	float invCount = 1.f / float(nu);
	for (int i = 0; i < nu+1; ++i)
	{
		cosTableU[i] = cosf(i * invCount * 2.f * PI);
		sinTableU[i] = sinf(i * invCount * 2.f * PI);
	}

	cosTableV = new float[nv+1];
	sinTableV = new float[nv+1];
	invCount = 1.f / float(nv);
	for (int i = 0; i < nv+1; ++i)
	{
		cosTableV[i] = cosf(i * invCount * PI);
		sinTableV[i] = sinf(i * invCount * PI);
	}
}

const Vector3 DomeLight::sampleLight(const unsigned int threadID, const Vector3 &from, const Vector3 &normal, const float time, const Scene &scene, const Vector3 &rVec, float &outSpec, bool isSecondary) const
{
	Ray sampleRay(threadID);
	HitInfo sampleHit;
	Vector3 direction = 0, tmpResult = 0, imageSample = 0;
	float e1, e2, fu, fv, theta, phi, cosTheta, sinTheta, cosPhi, sinPhi, tmpSpec = 0, pdfs[2], pdf;
	bool cutOff = false;
	int samplesDone = 0;
	float samplesDoneRecip = 1.0f;
	int numSamples = isSecondary ? 1 : m_numSamples;

	do
	{
		e1 = Scene::getRand(threadID); e2 = Scene::getRand(threadID);

		fu = uDistrib->sample(e1, &pdfs[0]);
		int u = ((int)fu == uDistrib->count) ? (int)fu-1 : (int)fu;
		fv = vDistribs[u]->sample(e2, &pdfs[1]);

		//theta = fv * vDistribs[u]->invCount * PI;
		//phi   = fu * uDistrib->invCount * 2.f * PI;

		cosTheta = cosTableV[int(fv)], sinTheta = sinTableV[int(fv)];
		sinPhi = sinTableU[int(fu)], cosPhi = cosTableU[int(fu)];

		direction = Vector3(-sinTheta*cosPhi, -cosTheta, -sinTheta*sinPhi);
		if (dot(normal, direction) < 0.0f) continue;

		pdf = (pdfs[0] * pdfs[1]) / (_2_PI2 * sinTheta);

		Vector3 E       = 0;
		float attenuate = 1.0f;
		imageSample = m_lightMap->getLookupXYZ3(direction);

		sampleHit.t = MIRO_TMAX;
		if (m_fastShadows)
		{
			sampleRay.set(threadID, from, direction, time, 1.001f, 0, 0, IS_SHADOW_RAY);		// Create shadow ray
			if (scene.trace(threadID, sampleHit, sampleRay, epsilon))					// Quick method, returns any hit
			{
				attenuate = 0.0f;
			}
		}
		else																// Full method, accounts for transparency effects
		{
			float distanceTraversed = 0.0f;
			sampleRay.set(threadID, from, direction, time, 1.001f, 0, 0, IS_PRIMARY_RAY);		// Create primary ray
			while (distanceTraversed < MIRO_TMAX && attenuate > epsilon)
			{
				if (scene.trace(threadID, sampleHit, sampleRay, epsilon))				
				{
					Vector3 hitN; sampleHit.getInterpolatedNormal(hitN);
					float nDL = dot(hitN, -direction);
					if (nDL > 0.0)											// Only attenuate on incoming direction
					{
						attenuate *= sampleHit.obj->m_material->refractAmt();
					}
					Vector3 newP = Vector3(sampleRay.o[0], sampleRay.o[1], sampleRay.o[2]) + sampleHit.t * direction;
					sampleRay.set(threadID, newP, direction, time, 1.001f, 0, 0, IS_PRIMARY_RAY);
					distanceTraversed += sampleHit.t;
				}
				else
				{
					distanceTraversed = MIRO_TMAX;
				}
			}
		}

		E = m_Gain * imageSample / pdf;										// Light irradiance for this sample
		samplesDone++;
		samplesDoneRecip = 1.0f / (float)samplesDone;

		cutOff = (E * samplesDoneRecip).average() < m_noiseThreshold;		// Stop sampling if contribution is below the noise threshold

		tmpResult += E * attenuate;
		tmpSpec   += dot(rVec, direction) * attenuate;

	}  while (samplesDone < numSamples && !cutOff);

	outSpec = tmpSpec * samplesDoneRecip;
	return tmpResult * samplesDoneRecip;
}