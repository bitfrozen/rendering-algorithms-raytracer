#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include <vector>

#define GET_NUMCHILD(a) (a)>>1
#define PUT_NUMCHILD(a) (a)<<1

ALIGN_SSE class BVH_Node
{
public:
	BVH_Node() {};
	~BVH_Node();
	ALIGN_SSE AABB* bBox;
	union {
		BVH_Node* Children;
		Object** objs;
	};
	union {
		u_int numChildren;
		u_int isLeaf;
	};
#ifdef USE_TRI_PACKETS
	ALIGN_SSE struct TriCache4 {
		Object* tris[4];
		union {float     Ax[4]; __m128     Ax4; };
		union {float     Ay[4]; __m128     Ay4; };
		union {float     Az[4]; __m128     Az4; };
		union {float edge0x[4]; __m128 edge0x4; };
		union {float edge0y[4]; __m128 edge0y4; };
		union {float edge0z[4]; __m128 edge0z4; };
		union {float edge1x[4]; __m128 edge1x4; };
		union {float edge1y[4]; __m128 edge1y4; };
		union {float edge1z[4]; __m128 edge1z4; };
	};
	TriCache4* triCache;
	void buildTriBundles(TriCache4* cacheAlloc);
#endif
	static unsigned int nodeCount, leafCount, maxDepth;

	bool intersect(HitInfo& result, const Ray& ray, float tMin = epsilon) const;

	void buildSAH(Object** objs, u_int numObjs, float* leftArea, float* rightArea);
	void partitionSweepSAH(Object** objs, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea);

	void buildBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, float* leftArea, float* rightArea, int* binIds);
	void partitionSweepBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, int* binIds);

	float calcSAHCost(int leftNum, float leftArea, int rightNum, float rightArea);
};

/*#define FIRST_NODE_IS_VALID		0x01
#define SECOND_NODE_IS_VALID	0x02
#define THIRD_NODE_IS_VALID		0x04
#define FOURTH_NODE_IS_VALID	0x08
#define FIRST_NODE_IS_LEAF		0x10
#define SECOND_NODE_IS_LEAF		0x20
#define THIRD_NODE_IS_LEAF		0x40
#define FOURTH_NODE_IS_LEAF		0x80

ALIGN_64 class QBVH_Node
{
public:
	QBVH_Node() {};
	~QBVH_Node();
	union {float bbMinX[4]; __m128 bbMinX4; };
	union {float bbMinY[4]; __m128 bbMinY4; };
	union {float bbMinZ[4]; __m128 bbMinZ4; };
	union {float bbMaxX[4]; __m128 bbMaxX4; };
	union {float bbMaxY[4]; __m128 bbMaxY4; };
	union {float bbMaxZ[4]; __m128 bbMaxZ4; };
	union {
		QBVH_Node* Children[4];
		BVH_Node::TriCache4* triCaches[4];
	};
	u_int flags;

	void build(BVH_Node* node, BVH_Node::TriCache4* cacheAlloc);
	void buildTriBundle(BVH_Node* node, BVH_Node::TriCache4* cacheAlloc, BVH_Node::TriCache4* triCache, int nodeNum);
};*/

class BVH
{
public:
    void build(Objects* objs);
	bool intersect(HitInfo& result, const Ray& ray, float tMin = epsilon) const;
	static unsigned long int rayBoxIntersections;
protected:
    Objects* m_objects;
	BVH_Node m_baseNode;
	//QBVH_Node m_baseQNode;
};

#endif // CSE168_BVH_H_INCLUDED
