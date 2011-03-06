#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include <vector>

#define GET_NUMCHILD(a)		(a>>3)
#define PUT_NUMCHILD(a)		(a<<3)
#define SET_IS_LEAF			0x01
#define IS_LEAF(a)			((a & 0x01) == true)
#define SET_AXIS(a)			(a & 0x03)<<1
#define PART_AXIS(a)		((a & 0x06)>>1)

ALIGN_SSE class BVH_Node
{
public:
	BVH_Node() {
		bBox = NULL;
		Children = NULL;
		numChildren = 0;
#ifdef USE_TRI_PACKETS
		triCache = NULL;
#endif
	}
	~BVH_Node();
	ALIGN_SSE AABB* bBox;
	union {
		BVH_Node* Children;
		Object** objs;
	};
	union {
		u_int numChildren;
		u_int flags;
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
		bool checkOut[4];
		bool doubleTake;
	};
	TriCache4* triCache;
	void buildTriBundles(TriCache4* cacheAlloc, bool reset = false);
#endif
	static unsigned int nodeCount, leafCount, maxDepth;
	const bool intersect(HitInfo &result, const Ray& ray, const float tMin = epsilon) const;

	void buildSAH(Object** objs, u_int numObjs, float* leftArea, float* rightArea);
	void partitionSweepSAH(Object** objs, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, u_int &bestAxis);

	void buildBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, float* leftArea, float* rightArea, int* binIds);
	void partitionSweepBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, int* binIds, u_int &bestAxis);

	float calcSAHCost(int leftNum, float leftArea, int rightNum, float rightArea);
};

#ifdef USE_QBVH

#define Q_GETAXIS_TOP(a)		(a & 0x003)
#define Q_GETAXIS_LEFT(a)		(a & 0x00C)>>2
#define Q_GETAXIS_RIGHT(a)		(a & 0x030)>>4
#define Q_SETAXIS_TOP(a)		(a & 0x003)
#define Q_SETAXIS_LEFT(a)		(a & 0x003)<<2
#define Q_SETAXIS_RIGHT(a)		(a & 0x003)<<4
#define Q_ISTOPLEAF(a)		    (a & 0x040)
#define Q_SET_TOPLEAF(a)		(a & 0x001)<<7
#define Q_ISLEFTLEAF(a)			(a & 0x080)
#define Q_SET_LEFTLEAF(a)		(a & 0x001)<<8
#define Q_ISRIGHTLEAF(a)		(a & 0x100)
#define Q_SET_RIGHTLEAF(a)		(a & 0x001)<<9
#define Q_ISNOLEAF(a)			(a & 0x200)
#define Q_SET_NOLEAF(a)			(a & 0x001)<<10

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
	bool flagsIsValid[4];
	bool flagsIsLeaf[4];
	u_int axisFlags;
	static unsigned int nodeCount, leafCount, maxDepth;

	const void intersect(HitInfo &result, const Ray& ray, __m128 &hit, const float tMin = epsilon) const;
	void build(BVH_Node* node, BVH_Node::TriCache4* cacheAlloc, bool reset = false);
	BVH_Node::TriCache4* buildTriBundle(BVH_Node* node, BVH_Node::TriCache4* cacheAlloc, int nodeNum);
	AABB getAABB();
	void getAABB(AABB* out);
};

#endif

class BVH
{
public:
	BVH() {
		m_baseNode = (BVH_Node*)_aligned_malloc(sizeof(BVH_Node), 16);
#ifdef USE_QBVH
		m_baseQNode = (QBVH_Node*)_aligned_malloc(sizeof(QBVH_Node), 64);
#endif
	}
	~BVH() {
		_aligned_free(m_baseNode);
#ifdef USE_QBVH
		_aligned_free(m_baseQNode);
#endif
	}
    void build(Objects* objs);

	void getAABB(AABB* outBox) {
#ifdef USE_QBVH
		m_baseQNode->getAABB(outBox);
#else
		outBox = m_baseNode->bBox;
#endif
	}
	AABB getAABB() {
#ifdef USE_QBVH
		return m_baseQNode->getAABB();
#else
		return *m_baseNode->bBox;
#endif
	}

	const bool intersect(const unsigned int threadID, HitInfo &result, const Ray& ray, const float tMin = epsilon) const;
	static unsigned long int rayBoxIntersections[256];
protected:
    Objects* m_objects;
	BVH_Node* m_baseNode;
#ifdef USE_QBVH
	ALIGN_64 mutable QBVH_Node* m_baseQNode;
#endif
};

#endif // CSE168_BVH_H_INCLUDED