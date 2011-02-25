#include "RectangleLight.h"


const Vector3 RectangleLight::getPointOnLight() const {
	float e1 = Scene::rands[Scene::randsIdx++];
	float e2 = Scene::rands[Scene::randsIdx++];
	e2 = (e2 > 0.99) ? 0.99 : e2;
	return m_v1 + e1*(m_v2-m_v1) + e2*(m_v3-m_v1);
}