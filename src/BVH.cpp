#include "BVH.h"
#include "Ray.h"
#include "Console.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include <xmmintrin.h>
#include <smmintrin.h>
#include <algorithm>
#include <time.h>

BVH_Node::~BVH_Node()
{
	delete bBox;
}

void BVH::build(Objects* objs)
{
	if (!use_BVH)
	{
		m_objects = objs;
		return;
	}
	u_int numObjs = objs->size();

	Object** BVHObjs = new Object*[numObjs];		// I use a double pointer to conserve memory.
	AABB* preCalcAABB = new AABB[numObjs];			// Pre-calculate all the AABBs to speed up BVH build.
	Vector3* centroids = new Vector3[numObjs];		// Pre-calculate all centroids
	for (int i = 0; i < numObjs; i++)				// This is used for the whole build, everything is done in place.
	{
		BVHObjs[i] = (*objs)[i];
		BVHObjs[i]->getAABB(&preCalcAABB[i]);
		centroids[i] = preCalcAABB[i].getCentroid();
	}
	float* leftArea = new float[num_bins];			// Scratch memory for the build
	float* rightArea = new float[num_bins];			// This way we prevent alloc/dealloc during build
	int* binIds = new int[numObjs];
	
	clock_t start = clock();
	m_baseNode.build(BVHObjs, preCalcAABB, centroids, numObjs, leftArea, rightArea, binIds);
	clock_t end = clock();

	delete[] leftArea;								// Clean up after ourselves...
	delete[] rightArea;
	delete[] binIds;
	delete[] preCalcAABB;
	delete[] centroids;
	
	printf("BVH Build time: %.4fs...\n", (end-start)/1000.f);
}

void BVH_Node::build(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, float* leftArea, float* rightArea, int* binIds)
{
	Vector3 nodeMin = Vector3(MIRO_TMAX);
	Vector3 nodeMax = Vector3(-MIRO_TMAX);
#ifndef NO_SSE
	bBox = (AABB*)_aligned_malloc(sizeof(AABB), 16);
#else
	bBox = new AABB;
#endif

	for (u_int i = 0; i < numObjs; i++)						// Get the AABB for this node
	{
		nodeMin.x = std::min(nodeMin.x, preCalcAABB[i].bbMin.x);
		nodeMin.y = std::min(nodeMin.y, preCalcAABB[i].bbMin.y);
		nodeMin.z = std::min(nodeMin.z, preCalcAABB[i].bbMin.z);
		nodeMax.x = std::max(nodeMax.x, preCalcAABB[i].bbMax.x);
		nodeMax.y = std::max(nodeMax.y, preCalcAABB[i].bbMax.y);
		nodeMax.z = std::max(nodeMax.z, preCalcAABB[i].bbMax.z);
	}
	*bBox = AABB(nodeMin, nodeMax);

	u_int partPt = 0;
	if (GET_NUMCHILD(numChildren = numObjs<<1) <= MAX_LEAF_SIZE)		// Make a leaf...
	{
		isLeaf |= true;
		this->objs = objs;
		qsort(objs, numObjs, sizeof(Object*), Object::sortByArea);		// Sort by area. This should speed up traversal a bit.
	}
	else
	{
		partitionSweep(objs, preCalcAABB, centroids, numObjs, partPt, leftArea, rightArea, binIds);		// Find the best place to split the node
		numChildren = NODE_SIZE<<1;										// Check out BVH.h, the first of numChildren is used by isLeaf.
		isLeaf |= false;

		Children = new BVH_Node[2];
		Object** leftObjs = objs;				// Get a new pointer so we don't have to do so much pointer arithmetic inside partitionSweep
		AABB* leftAABB = preCalcAABB;			//
		Vector3* leftCentroids = centroids;		//
		u_int leftNum = partPt+1;				// [0 .. partPt]
		Object** rightObjs = objs+(partPt+1);	// Same as before
		AABB* rightAABB = preCalcAABB+(partPt+1);
		Vector3* rightCentroids = centroids+(partPt+1);
		u_int rightNum = numObjs-partPt-1;		// [partPt+1 .. numObjs-1]

		Children[0].build(leftObjs, leftAABB, leftCentroids, leftNum, leftArea, rightArea, binIds);	// Build left and right child nodes.
		Children[1].build(rightObjs, rightAABB, rightCentroids, rightNum, leftArea, rightArea, binIds);
		return;
	}
}

void BVH_Node::partitionSweep(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, int* binIds)
{
	float bestCost  = MIRO_TMAX;
	float thisCost  = 0;
	int bestAxis	= -1;
	int binPart	= 0;

	if (use_Bins)		// Do a binned build. Loosely following some guidelines from Wald's 2007 paper http://www.sci.utah.edu/~wald/Publications/2007/fastBuild/download/fastbuild.pdf
	{									
		AABB binBounds = AABB();

		for (int i = 0; i < numObjs; i++)		// Find the bounds defined by the objects' centroids
		{
			binBounds.grow(centroids[i]);
		}
		float length[3] = {binBounds.bbMax.x - binBounds.bbMin.x,		// Get dimensions of the bounds.
						   binBounds.bbMax.y - binBounds.bbMin.y,		// We only test along the longest axis
						   binBounds.bbMax.z - binBounds.bbMin.z};
		int axis[3];
		if (length[0] > length[1])										// Check for longest axis.
		{																// We store the order of lengths in case we
			if (length[0] > length[2])									// want to test all axis' at some point.
			{
				axis[0] = 0;
				axis[1] = (length[1] > length[2]) ? 1 : 2;
				axis[2] = (length[1] > length[2]) ? 2 : 1;
			}
			else
			{
				axis[0] = 2;
				axis[1] = 0;
				axis[2] = 1;
			}
		}
		else
		{
			if (length[1] > length[2])
			{
				axis[0] = 1;
				axis[1] = (length[0] > length[2]) ? 0 : 2;
				axis[2] = (length[0] > length[2]) ? 2 : 0;
			}
			else
			{
				axis[0] = 2;
				axis[1] = 1;
				axis[2] = 0;
			}
		}

		float kl = (float)num_bins*(1.0f-epsilon) / length[axis[0]];		// Coefficients for binning the objects
		float ko = binBounds.bbMin[axis[0]];								// Bin(centroid[i]) = num_bins*(centroid[i][axis] - bounds.min[axis]) / length[axis]

		AABB binBBs[num_bins];
		int numTris[num_bins];
		for (int i = 0; i < num_bins; i++)
		{
			numTris[i] = 0;
		}

		float newCentroid;
		for (int i = 0; i < numObjs; i++)		// Bin the objects. Store the number of objects in each bin.
		{
			newCentroid = centroids[i][axis[0]];
			binIds[i] = kl*(newCentroid - ko);

			binBBs[binIds[i]] = AABB(binBBs[binIds[i]], preCalcAABB[i]);
			numTris[binIds[i]]++;
		}

		// sweep from left
		AABB tempBBox;
		for (int i = 0; i < num_bins-1; i++)
		{
			tempBBox = AABB(tempBBox, binBBs[i]);
			leftArea[i] = tempBBox.getArea();		// Surface area of AABB for the left i bins.
		}

		// sweep from right
		tempBBox = AABB();
		int tempNum = 0;
		for (int i = num_bins-1; i > 0; i--)	
		{
			tempNum += numTris[i];
			tempBBox = AABB(tempBBox, binBBs[i]);
			rightArea[i] = tempBBox.getArea();											// Surface area of AABB for right bins.
			thisCost = (numObjs-tempNum)*leftArea[i-1] + tempNum*rightArea[i];		// Calculate the cost of this partition.
			if (thisCost < bestCost)												// If this is better than before, use this partition.
			{
				bestCost = thisCost;
				binPart = i;
			}
		}
		int revIdx = numObjs-1;				
		for (int i = 0; i < numObjs; i++)		// Loose sorting. We just ensure that all objects from the left bins are first in the objs pointer.
		{
			if (binIds[i] >= binPart)
			{
				while (binIds[revIdx] >= binPart) {revIdx--;}
				if (revIdx <= i)
				{
					partPt = i-1;
					return;
				}
				AABB tmpBox = preCalcAABB[i];
				preCalcAABB[i] = preCalcAABB[revIdx];
				preCalcAABB[revIdx] = tmpBox;
				Vector3 tmpVec = centroids[i];
				centroids[i] = centroids[revIdx];
				centroids[revIdx] = tmpVec;
				Object* tmp = objs[i];
				objs[i] = objs[revIdx];
				objs[revIdx--] = tmp;
			}
		}
	}
	else		// If there's fewer objects than bins, we just do full sweeps. 
	{			// Using ideas from Wald's "Ray Tracing Deformable Scenes Using Dynamic Bounding Volume Hierarchies" (Siggraph, 2007).
		
		// sweep in the x direction. First we sort (note we only sort the range we need).
		qsort(objs, numObjs, sizeof(Object*), Object::sortByXComponent);

		// sweep from left
		AABB tempBBox;
		AABB newBBox;

		leftArea[0] = MIRO_TMAX;
		for (int i = 1; i < numObjs; i++)
		{
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();		// Surface area of AABB for left side of node.
		}

		// sweep from right
		tempBBox = AABB();
		rightArea[numObjs-1] = MIRO_TMAX;
		for (int i = numObjs-2; i >= 0; i--)
		{
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();									// Surface area of AABB for right side of node.
			thisCost = (i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i];		// Calculate the cost of this partition.
			if (thisCost < bestCost)										// If this is better than before, use this partition.
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
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();
		}

		// sweep from right
		tempBBox = AABB();
		rightArea[numObjs-1] = MIRO_TMAX;
		for (int i = numObjs-2; i >= 0; i--)
		{
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();
			thisCost = (i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i];
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
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();
		}

		// sweep from right
		tempBBox = AABB();
		rightArea[numObjs-1] = MIRO_TMAX;
		for (int i = numObjs-2; i >= 0; i--)
		{
			objs[i]->getAABB(&newBBox);
			tempBBox = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();
			thisCost = (i+1)*leftArea[i] + (numObjs-i-1)*rightArea[i];
			if (thisCost < bestCost)
			{
				bestCost = thisCost;
				partPt = i;
				bestAxis = 3;
			}
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
	}
}

int Object::sortByArea(const void* s1, const void* s2)
{
	float left, right;
	left = (*(Object**)s1)->getAABB().getArea();
	right = (*(Object**)s2)->getAABB().getArea();

	if (left < right)
	{
		return -1;
	}
	else if (left > right)
	{
		return 1;
	}
	return 0;
}

int Object::sortByXComponent(const void* s1, const void* s2)
{
	Vector3 left, right;
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

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
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

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
	left = (*(Object**)s1)->getAABB().getCentroid();
	right = (*(Object**)s2)->getAABB().getCentroid();

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

bool BVH_Node::intersect(HitInfo& result, const Ray& ray, float tMin /* = epsilon */, float tMax /* = MIRO_TMAX */) const
{
	bool hit = false;

#ifndef NO_SSE														// Implementation of Manny Ko's algorithm, from http://tog.acm.org/resources/RTNews/html/rtnv23n1.html#art7
																	// Had to guess as to a few details, but works great.
	const __m128 rayorig = _mm_load_ps(&ray.ox);					
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
