#include "BVH.h"
#include "Ray.h"
#include "Console.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include <xmmintrin.h>
#include <smmintrin.h>
#include <algorithm>

void
BVH::build(Objects* objs)
{
	if (!use_BVH)
	{
		m_objects = objs;
		return;
	}
	u_int numObjs = objs->size();

	Object** BVHObjs = new Object*[numObjs];		// I use a double pointer to conserve memory.
	for (int i = 0; i < numObjs; i++)				// This is used for the whole build, everything is done in place.
	{
		BVHObjs[i] = (*objs)[i];
	}
	m_baseNode.build(BVHObjs, numObjs);

	// Clean up BVH Precalc memory used by meshes
	Objects::const_iterator objIter;
	Triangle* triMesh;
	for (objIter = objs->begin(); objIter != objs->end(); objIter++)
	{
		if (triMesh = dynamic_cast<Triangle*>(*objIter))
		{
			triMesh->cleanBVHMem();
		}
	}
}

bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
	if (!use_BVH)
	{
		bool hit = false;

		for (size_t i = 0; i < m_objects->size(); ++i)
		{
			if ((*m_objects)[i]->intersect(minHit, ray, tMin, tMax))
			{
				hit = true;
			}
		}    
		return hit;
	}
	return m_baseNode.intersect(minHit, ray, tMin, tMax);
}

void BVH_Node::build(Object** objs, u_int numObjs)//, centroidMap* partitionMap)
{
	Vector3 nodeMin = Vector3(MIRO_TMAX);
	Vector3 nodeMax = Vector3(-MIRO_TMAX);	
	AABB tempBox;
	bBox = (AABB*)_aligned_malloc(sizeof(AABB), 16);

	for (u_int i = 0; i < numObjs; i++)						// Get the AABB for this node
	{
		objs[i]->getAABB(tempBox);
		nodeMin.x = std::min(nodeMin.x, tempBox.bbMin.x);
		nodeMin.y = std::min(nodeMin.y, tempBox.bbMin.y);
		nodeMin.z = std::min(nodeMin.z, tempBox.bbMin.z);
		nodeMax.x = std::max(nodeMax.x, tempBox.bbMax.x);
		nodeMax.y = std::max(nodeMax.y, tempBox.bbMax.y);
		nodeMax.z = std::max(nodeMax.z, tempBox.bbMax.z);
	}
	*bBox = AABB(nodeMin, nodeMax);

	u_int partPt = 0;
	if (GET_NUMCHILD(numChildren = numObjs<<1) <= MAX_LEAF_SIZE)		// Make a leaf...
	{
		isLeaf |= true;
		this->objs = objs;
	}
	else
	{
		partitionSweep(objs, numObjs, partPt, 1.0f/bBox->getArea());	// Find the best place to split the node
		numChildren = NODE_SIZE<<1;										// Check out BVH.h, the first of numChildren is used by isLeaf.
		isLeaf |= false;

		Children = new BVH_Node[2];
		Object** leftObjs = objs;				// Get a new pointer so we don't have to do so much pointer arithmetic inside partitionSweep
		u_int leftNum = partPt+1;				// [0 .. partPt]
		Object** rightObjs = objs+(partPt+1);	// Same as before
		u_int rightNum = numObjs-partPt-1;		// [partPt+1 .. numObjs-1]

		Children[0].build(leftObjs, leftNum);	// Build left and right child nodes.
		Children[1].build(rightObjs, rightNum);
		return;
	}
}

bool BVH_Node::partitionSweep(Object** objs, u_int numObjs, u_int& partPt, float bbSAInv)
{
	float bestCost  = MIRO_TMAX;
	float thisCost  = 0;
	int bestAxis	= -1;
	float* leftArea = new float[numObjs];
	float* rightArea = new float[numObjs];

	// sweep in the x direction. First we sort (note we only sort the range we need).
	qsort(objs, numObjs, sizeof(Object*), Object::sortByXComponent);

	// sweep from left
	AABB tempBBox = AABB();
	AABB newBBox;
	leftArea[0] = MIRO_TMAX;
	for (int i = 1; i < numObjs; i++)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		leftArea[i] = tempBBox.getArea();		// Surface area of AABB for left side of node.
	}

	// sweep from right
	tempBBox = AABB();
	rightArea[numObjs-1] = tempBBox.getArea();
	for (int i = numObjs-2; i >= 0; i--)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		rightArea[i] = tempBBox.getArea();														// Surface area of AABB for right side of node.
		thisCost = 2*Tbox + bbSAInv*Ttri*((i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i]);		// Calculate the cost of this partition.
		if (thisCost < bestCost)																// If this is better than before, use this partition.
		{
			bestCost = thisCost;
			partPt = i;
			bestAxis = 1;
		}
	}

	// sweep in the y direction
	qsort(objs, numObjs, sizeof(Object*), Object::sortByYComponent);

	// sweep from left
	tempBBox = AABB();
	leftArea[0] = MIRO_TMAX;
	for (int i = 1; i < numObjs; i++)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		leftArea[i] = tempBBox.getArea();
	}

	// sweep from right
	tempBBox = AABB();
	rightArea[numObjs-1] = tempBBox.getArea();
	for (int i = numObjs-2; i >= 0; i--)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		rightArea[i] = tempBBox.getArea();
		thisCost = 2*Tbox + bbSAInv*Ttri*((i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i]);
		if (thisCost < bestCost)
		{
			bestCost = thisCost;
			partPt = i;
			bestAxis = 2;
		}
	}

	// sweep in the z direction
	qsort(objs, numObjs, sizeof(Object*), Object::sortByZComponent);

	// sweep from left
	tempBBox = AABB();
	leftArea[0] = MIRO_TMAX;
	for (int i = 1; i < numObjs; i++)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		leftArea[i] = tempBBox.getArea();
	}

	// sweep from right
	tempBBox = AABB();
	rightArea[numObjs-1] = tempBBox.getArea();
	for (int i = numObjs-2; i >= 0; i--)
	{
		objs[i]->getAABB(newBBox);
		tempBBox = AABB(tempBBox, newBBox);
		rightArea[i] = tempBBox.getArea();
		thisCost = 2*Tbox + bbSAInv*Ttri*((i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i]);
		if (thisCost < bestCost)
		{
			bestCost = thisCost;
			partPt = i;
			bestAxis = 3;
		}
	}
	delete[] leftArea;	// Clean up...
	delete[] rightArea;

	if (bestAxis == -1) 
	{
		return true; // make leaf node, no better partition found
	}

	switch (bestAxis)
	{
	case 1 :
		qsort(objs, numObjs, sizeof(Object*), Object::sortByXComponent);
		break;
	case 2 :
		qsort(objs, numObjs, sizeof(Object*), Object::sortByYComponent);
		break;
	}
	return false; // make child nodes
}

int Object::sortByXComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	(*(Object**)s1)->getCentroid(left);
	(*(Object**)s2)->getCentroid(right);

	if (left.x < right.x)
	{
		return -1;
	}
	else if (left.x > right.x)
	{
		return 1;
	}
	return 0;
}

int Object::sortByYComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	(*(Object**)s1)->getCentroid(left);
	(*(Object**)s2)->getCentroid(right);

	if (left.y < right.y)
	{
		return -1;
	}
	else if (left.y > right.y)
	{
		return 1;
	}
	return 0;
}

int Object::sortByZComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	(*(Object**)s1)->getCentroid(left);
	(*(Object**)s2)->getCentroid(right);

	if (left.z < right.z)
	{
		return -1;
	}
	else if (left.z > right.z)
	{
		return 1;
	}
	return 0;
}

bool BVH_Node::intersect(HitInfo& result, const Ray& ray, float tMin /* = epsilon */, float tMax /* = MIRO_TMAX */) const
{
	bool hit = false;

#ifndef NO_SSE														// Implementation of Manny Ko's algorithm, from http://tog.acm.org/resources/RTNews/html/rtnv23n1.html#art7
	const __m128 rayorig = _mm_load_ps(&ray.ox);					// Had to guess as to a few details, but works great.
	const __m128 rayinvdir = _mm_load_ps(&ray.idx);
	const __m128 minPoint = _mm_load_ps(&bBox->bbMin.x);
	const __m128 maxPoint = _mm_load_ps(&bBox->bbMax.x);
	__m128 t0v, t1v;
	
	t0v = _mm_mul_ps(_mm_sub_ps(minPoint, rayorig), rayinvdir);
	t1v = _mm_mul_ps(_mm_sub_ps(maxPoint, rayorig), rayinvdir);

	__m128 t0v_n, t1v_n;
	t0v_n = _mm_min_ps(t0v, t1v);
	t1v_n = _mm_max_ps(t0v, t1v);

	__m128 interval_min, interval_max, t0, t1, xmm0, xmm1;			// Horizontal max, can't find a single SSE instruction that will do this, though Manny Ko seems to imply there is one...
	xmm0 = t0v_n;
	xmm1 = _mm_movehl_ps( xmm0, xmm0 );
	xmm0 = _mm_max_ps( xmm0, xmm1 );
	xmm1 = _mm_shuffle_ps( xmm0, xmm0, _MM_SHUFFLE(1,1,1,1) );
	t0   = _mm_max_ss( xmm0, xmm1 );

	xmm0 = t1v_n;													// Horizontal min, same deal as before...
	xmm1 = _mm_shuffle_ps( xmm0, xmm0, _MM_SHUFFLE(2,2,2,2) );
	xmm0 = _mm_min_ps( xmm0, xmm1 );
	xmm1 = _mm_shuffle_ps( xmm0, xmm0, _MM_SHUFFLE(1,1,1,1) );
	t1   = _mm_min_ss( xmm0, xmm1 );

	interval_min = _mm_max_ss(t0, _mm_load_ss(&tMin)); // tMin
	interval_max = _mm_min_ss(t1, _mm_load_ss(&tMax)); // tMax

	if (_mm_comile_ss(interval_min, interval_max))
	{
		if (isLeaf & true)											// If this is a leaf, then we need to intersect the objects inside...
		{
			for (int i = 0; i < GET_NUMCHILD(numChildren); i++)
			{
				if (objs[i]->intersect(result, ray, tMin, tMax))
				{
					hit = true;
				}
			}
		}
		else														// Otherwise, recurse on the child nodes...
		{
			for (int i = 0; i < GET_NUMCHILD(numChildren); i++)
			{				
				if (Children[i].intersect(result, ray, tMin, tMax))
				{
					hit = true;
				}
			}
		}
	}
#else																// Standard slabs test, following Williams et al. http://portal.acm.org/citation.cfm?id=1198748
	float bbmin, bbmax, tymin, tymax, tzmin, tzmax;	

	bbmin = ((&bBox->bbMin)[ray.bounces_flags & SIGN_X_FLAG].x - ray.ox) * ray.idx;
	bbmax = ((&bBox->bbMin)[1-(ray.bounces_flags & SIGN_X_FLAG)].x - ray.ox) * ray.idx;
	tymin = ((&bBox->bbMin)[GET_Y_SIGN(ray.bounces_flags & SIGN_Y_FLAG)].y - ray.oy) * ray.idy;
	tymax = ((&bBox->bbMin)[1-(GET_Y_SIGN(ray.bounces_flags & SIGN_Y_FLAG))].y - ray.oy) * ray.idy;

	if ((bbmin > tymax) || (tymin > bbmax))
		return false;

	if (tymin > bbmin)
		bbmin = tymin;
	if (tymax < bbmax)
		bbmax = tymax;

	tzmin = ((&bBox->bbMin)[GET_Z_SIGN(ray.bounces_flags & SIGN_Z_FLAG)].z - ray.oz) * ray.idz;
	tzmax = ((&bBox->bbMin)[1-(GET_Z_SIGN(ray.bounces_flags & SIGN_Z_FLAG))].z - ray.oz) * ray.idz;

	if ((bbmin > tzmax) || (tzmin > bbmax))
		return false;
	if (tzmin > bbmin)
		bbmin = tzmin;
	if (tzmax < bbmax)
		bbmax = tzmax;

	if ((bbmin <= tMax) && (bbmax >= tMin));
	{
		if (isLeaf & true)
		{
			for (int i = 0; i < GET_NUMCHILD(numChildren); i++)
			{
				if (objs[i]->intersect(result, ray, tMin, tMax))
				{
					hit = true;
				}
			}
		}
		else
		{
			for (int i = 0; i < GET_NUMCHILD(numChildren); i++)
			{				
				if (Children[i].intersect(result, ray, tMin, tMax))
				{
					hit = true;
				}
			}
		}
	}
#endif
	return hit;
}