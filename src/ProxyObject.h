#ifndef CSE168_PROXYOBJECT_H_INCLUDED
#define CSE168_PROXYOBJECT_H_INCLUDED

#include "Object.h"
#include "BVH.h"
#include "Material.h"
#include "ProxyMatrix.h"

using namespace std;

ALIGN_SSE class ProxyObject : public Object
{
public:
	ProxyObject(Objects* m = NULL, BVH* b = NULL, Matrix4x4 &t = Matrix4x4());
	~ProxyObject() {};

	void setDisplayNum(int n) {m_displayNum = n;}
	void renderGL();
	void preCalc();

	const bool intersect(const unsigned int threadID, HitInfo &result, const Ray& ray, const float tMin = epsilon);
	const ProxyMatrix& getMatrix() const {return *m_Matrix;}

	void getAABB(AABB* outBox);
	AABB getAABB();

	static void setupProxy(TriangleMesh* mesh = NULL, Material* mat = NULL, Objects* m = NULL, BVH* b = NULL);
	static void setupMultiProxy(TriangleMesh* mesh[], int numObjs, Material* mat[], Objects* m, BVH* b);

private:
	Objects* m_objects;
	BVH* m_BVH;
	ProxyMatrix* m_Matrix;
	int m_displayNum;
};

#endif // CSE168_PROXYOBJECT_H_INCLUDED