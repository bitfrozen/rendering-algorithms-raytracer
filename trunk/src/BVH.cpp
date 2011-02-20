#include "BVH.h"
#include "Ray.h"
#include "Console.h"
#include "TriangleMesh.h"
#include "SSE.h"
#include <algorithm>
#include <time.h>
#include <omp.h>
#include "Object.h"

using namespace std;

unsigned int BVH_Node::nodeCount = 0;
unsigned int BVH_Node::leafCount = 0;
unsigned int BVH_Node::maxDepth  = 0;

BVH_Node::~BVH_Node()
{
	delete bBox;
#ifdef USE_TRI_PACKETS
	delete triCache;
#endif
}

void BVH::build(Objects* objs)
{
	if (!use_BVH)
	{
		m_objects = objs;
		return;
	}

	u_int numObjs = objs->size();

#ifdef USE_BINS

	Object**   BVHObjs = new Object*[numObjs];		// Use a double pointer to alias into scene object vector container.
	float*    leftArea = new   float[128+NUM_BINS];	// Scratch memory for the build
	float*   rightArea = new   float[128+NUM_BINS];	// This way we prevent alloc/dealloc during build
	Vector3* centroids = new Vector3[numObjs];		// Pre-calculate all centroids
	AABB*  preCalcAABB = new    AABB[numObjs];		// Pre-calculate all the AABBs to speed up BVH build.
	int*        binIds = new     int[numObjs];		// Temporary memory to hold the binIDs of the objects.
	int*	   numTris = new     int[NUM_BINS];		// Scratch memory for the build
	AABB*	    binBBs = new	AABB[NUM_BINS];		// Scratch memory for the build

	for (int i = 0; i < numObjs; i++)				// This is used for the whole build, everything is done in place.
	{
		BVHObjs[i] = (*objs)[i];
		centroids[i] = BVHObjs[i]->getAABB().getCentroid();
		BVHObjs[i]->getAABB(&preCalcAABB[i]);
	}

	clock_t start = clock();

	m_baseNode.buildBin(BVHObjs, preCalcAABB, centroids, numObjs, leftArea, rightArea, binIds);

#ifdef USE_TRI_PACKETS
	BVH_Node::TriCache4* cacheAlloc = (BVH_Node::TriCache4*)_aligned_malloc(sizeof(BVH_Node::TriCache4)*BVH_Node::leafCount, 16);
	m_baseNode.buildTriBundles(cacheAlloc);
#endif

	clock_t end = clock();

	delete[] binBBs;
	delete[] numTris;
	delete[] binIds;								// Clean up after ourselves...
	delete[] preCalcAABB;
	delete[] centroids;
	delete[] rightArea;
	delete[] leftArea;

#else

	Object**   BVHObjs = new Object*[numObjs];		// Use a double pointer to alias into scene object vector container.
	float*    leftArea = new   float[numObjs];		// Scratch memory for the build
	float*   rightArea = new   float[numObjs];		// This way we prevent alloc/dealloc during build

	for (int i = 0; i < numObjs; i++)				// This is used for the whole build, everything is done in place.
	{
		BVHObjs[i] = (*objs)[i];
	}

	clock_t start = clock();

	m_baseNode.buildSAH(BVHObjs, numObjs, leftArea, rightArea);

#ifdef USE_TRI_PACKETS
	BVH_Node::TriCache4* cacheAlloc = (BVH_Node::TriCache4*)_aligned_malloc(sizeof(BVH_Node::TriCache4)*BVH_Node::leafCount, 16);
	m_baseNode.buildTriBundles(cacheAlloc);
#endif

	clock_t end = clock();

	delete[] rightArea;
	delete[] leftArea;

#endif

	printf("BVH Build time: %.4fs...\n",	(end-start)/1000.f);
	printf("Number of nodes: %d...\n",		BVH_Node::nodeCount);
	printf("Number of leaves: %d...\n",		BVH_Node::leafCount);
	printf("Maximum depth: %d...\n",		BVH_Node::maxDepth);
	printf("Average faces/leaf: %.4f...\n", numObjs/(float)BVH_Node::leafCount);
}

#ifdef USE_TRI_PACKETS
void BVH_Node::buildTriBundles(TriCache4* cacheAlloc)
{
	static int nodeNum = 0;
	if (!(isLeaf & true))
	{
		Children[0].buildTriBundles(cacheAlloc);
		Children[1].buildTriBundles(cacheAlloc);
	}
	else
	{
		triCache = &cacheAlloc[nodeNum++];
		for (int i = 0; i < 4; i++)
		{
			triCache->Ax[i] = triCache->Ay[i] = triCache->Az[i] \
				= triCache->edge0x[i] = triCache->edge0y[i] = triCache->edge0z[i] \
				= triCache->edge1x[i] = triCache->edge1y[i] = triCache->edge1z[i] = 0;
			triCache->tris[i] = 0;
		}
		for (int i = 0; i < GET_NUMCHILD(numChildren); i++)
		{
			TriangleMesh::TupleI3 ti3;
			ti3                 = objs[i]->m_mesh->m_vertexIndices[objs[i]->m_index];
			triCache->Ax[i]     = objs[i]->m_mesh->m_vertices[ti3.x].x;
			triCache->Ay[i]     = objs[i]->m_mesh->m_vertices[ti3.x].y;
			triCache->Az[i]     = objs[i]->m_mesh->m_vertices[ti3.x].z;
			triCache->edge0x[i] = objs[i]->m_mesh->m_vertices[ti3.y].x - objs[i]->m_mesh->m_vertices[ti3.x].x;
			triCache->edge0y[i] = objs[i]->m_mesh->m_vertices[ti3.y].y - objs[i]->m_mesh->m_vertices[ti3.x].y;
			triCache->edge0z[i] = objs[i]->m_mesh->m_vertices[ti3.y].z - objs[i]->m_mesh->m_vertices[ti3.x].z;
			triCache->edge1x[i] = objs[i]->m_mesh->m_vertices[ti3.z].x - objs[i]->m_mesh->m_vertices[ti3.x].x;
			triCache->edge1y[i] = objs[i]->m_mesh->m_vertices[ti3.z].y - objs[i]->m_mesh->m_vertices[ti3.x].y;
			triCache->edge1z[i] = objs[i]->m_mesh->m_vertices[ti3.z].z - objs[i]->m_mesh->m_vertices[ti3.x].z;
			triCache->tris[i]   = objs[i];
		}
	}
}
#endif

void BVH_Node::buildBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, float* leftArea, float* rightArea, int* binIds)
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
		nodeMin.x = min(nodeMin.x, objs[i]->getAABB().bbMin.x);
		nodeMin.y = min(nodeMin.y, objs[i]->getAABB().bbMin.y);
		nodeMin.z = min(nodeMin.z, objs[i]->getAABB().bbMin.z);
		nodeMax.x = max(nodeMax.x, objs[i]->getAABB().bbMax.x);
		nodeMax.y = max(nodeMax.y, objs[i]->getAABB().bbMax.y);
		nodeMax.z = max(nodeMax.z, objs[i]->getAABB().bbMax.z);
	}
	*bBox = AABB(nodeMin, nodeMax);

	u_int partPt = 0;
	if (numObjs <= MAX_LEAF_SIZE)		// Make a leaf...
	{
		numChildren = PUT_NUMCHILD(numObjs);
		isLeaf |= true;
		nodeCount++;
		leafCount++;

		this->objs = objs;
		qsort(objs, numObjs, sizeof(Object*), Object::sortByArea);
	}
	else
	{
		static int curDepth = 0;

		curDepth++;
		if (curDepth > maxDepth) maxDepth = curDepth;
		nodeCount += 2;

		partitionSweepBin(objs, preCalcAABB, centroids, numObjs, partPt, leftArea, rightArea, binIds);		// Find the best place to split the node

		numChildren = PUT_NUMCHILD(NODE_SIZE);															// Check out BVH.h, the first bits of numChildren is used by isLeaf.
		isLeaf |= false;

		Children = (BVH_Node*)_aligned_malloc(sizeof(BVH_Node)*2, 16);

		Object** leftObjs       = objs;						// Get a new pointer so we don't have to do so much pointer arithmetic inside partitionSweep
		AABB* leftAABB          = preCalcAABB;				//
		Vector3* leftCentroids  = centroids;				//
		u_int leftNum           = partPt+1;					// [0 .. partPt]

		Object** rightObjs      = objs+(partPt+1);			// Same as before
		AABB* rightAABB         = preCalcAABB+(partPt+1);	//
		Vector3* rightCentroids = centroids+(partPt+1);		//
		u_int rightNum          = numObjs-partPt-1;			// [partPt+1 .. numObjs-1]

		Children[0].buildBin(leftObjs,   leftAABB,  leftCentroids,  leftNum, leftArea, rightArea, binIds);	// Build left and right child nodes.
		Children[1].buildBin(rightObjs, rightAABB, rightCentroids, rightNum, leftArea, rightArea, binIds);
		curDepth--;
		return;
	}
}

void BVH_Node::partitionSweepBin(Object** objs, AABB* preCalcAABB, Vector3* centroids, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea, int* binIds)
{
	float bestCost  = INFINITY;
	float thisCost  = 0;
	int bestAxis	= -1;
	int binPart		= 0;

#ifdef USE_TRI_PACKETS
	if (numObjs >= 128)
#endif
	{
		// Do a binned build. Loosely following some guidelines from Wald's 2007 paper http://www.sci.utah.edu/~wald/Publications/2007/fastBuild/download/fastbuild.pdf								
		AABB binBounds = AABB();

		for (int i = 0; i < numObjs; i++)		// Find the bounds defined by the objects' centroids
		{
			binBounds.grow(centroids[i]);
		}
		float length[3] = {binBounds.bbMax.x - binBounds.bbMin.x,		// Get dimensions of the bounds.
			binBounds.bbMax.y - binBounds.bbMin.y,		// We test along all axis to find the best partition
			binBounds.bbMax.z - binBounds.bbMin.z};

		for (int axis = 0; axis < 3; axis++)
		{
			if (length[axis] < epsilon)
			{
				break;
			}
			float kl = (float)NUM_BINS*(1.0f-epsilon) / length[axis];		// Coefficients for binning the objects
			float ko = binBounds.bbMin[axis];								// Bin(centroid[i]) = num_bins*(centroid[i][axis] - bounds.min[axis]) / length[axis]

			AABB binBBs[NUM_BINS];
			int numTris[NUM_BINS];
			for (int i = 0; i < NUM_BINS; i++)
			{
				numTris[i] = 0;
			}

			float newCentroid;
			for (int i = 0; i < numObjs; i++)		// Bin the objects. Store the number of objects in each bin.
			{
				newCentroid			= centroids[i][axis];
				binIds[i]			= kl*(newCentroid - ko);

				binBBs[binIds[i]]	= AABB(binBBs[binIds[i]], preCalcAABB[i]);
				numTris[binIds[i]]++;
			}

			// sweep from left
			AABB tempBBox;
			for (int i = 0; i < NUM_BINS-1; i++)
			{
				tempBBox	= AABB(tempBBox, binBBs[i]);
				leftArea[i] = tempBBox.getArea();		// Surface area of AABB for the left i bins.
			}

			// sweep from right
			tempBBox = AABB();
			int tempNum = 0;
			for (int i = NUM_BINS-1; i > 0; i--)	
			{
				tempNum		+= numTris[i];
				tempBBox	 = AABB(tempBBox, binBBs[i]);
				rightArea[i] = tempBBox.getArea();														// Surface area of AABB for right bins.
				thisCost	 = calcSAHCost( (numObjs-tempNum), leftArea[i-1], tempNum, rightArea[i]);	// Calculate the cost of this partition.
				if (thisCost < bestCost)																// If this is better than before, use this partition.
				{
					bestCost = thisCost;
					binPart = i;
					bestAxis = axis;
				}
			}
		}
		float kl = (float)NUM_BINS*(1.0f-epsilon) / length[bestAxis];		// Coefficients for binning the objects
		float ko = binBounds.bbMin[bestAxis];								// Bin(centroid[i]) = num_bins*(centroid[i][axis] - bounds.min[axis]) / length[axis]

		float newCentroid;
		for (int i = 0; i < numObjs; i++)		// Bin the objects (we need to do this again to sort below)
		{
			newCentroid	= centroids[i][bestAxis];
			binIds[i]   = kl*(newCentroid - ko);
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
#ifdef USE_TRI_PACKETS
	else
	{
		// Using ideas from Wald's "Ray Tracing Deformable Scenes Using Dynamic Bounding Volume Hierarchies" (Siggraph, 2007).	
		// sweep in the x direction. First we sort (note we only sort the range we need).
		qsort(objs, numObjs, sizeof(Object*), Object::sortByXComponent);

		// sweep from left
		AABB tempBBox;
		AABB newBBox;

		leftArea[0] = INFINITY;
		for (int i = 1; i < numObjs; i++)
		{
			newBBox     = objs[i]->getAABB();
			tempBBox    = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();		// Surface area of AABB for left side of node.
		}

		// sweep from right
		tempBBox             = AABB();
		rightArea[numObjs-1] = INFINITY;
		for (int i = numObjs-2; i >= 0; i--)
		{
			newBBox      = objs[i]->getAABB();
			tempBBox     = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();												// Surface area of AABB for right side of node.
			thisCost     = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);	// Calculate the cost of this partition.
			if (thisCost < bestCost)														// If this is better than before, use this partition.
			{
				bestCost = thisCost;
				partPt   = i;
				bestAxis = 0;
			}
		}

		// sweep in the y direction
		qsort(objs, numObjs, sizeof(Object*), Object::sortByYComponent);

		// sweep from left
		tempBBox    = AABB();
		leftArea[0] = INFINITY;
		for (int i = 1; i < numObjs; i++)
		{
			newBBox     = objs[i]->getAABB();
			tempBBox    = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();
		}

		// sweep from right
		tempBBox             = AABB();
		rightArea[numObjs-1] = INFINITY;
		for (int i = numObjs-2; i >= 0; i--)
		{
			newBBox      = objs[i]->getAABB();
			tempBBox     = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();
			thisCost     = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);
			if (thisCost < bestCost)
			{
				bestCost = thisCost;
				partPt   = i;
				bestAxis = 1;
			}
		}

		// sweep in the z direction
		qsort(objs, numObjs, sizeof(Object*), Object::sortByZComponent);

		// sweep from left
		tempBBox    = AABB();
		leftArea[0] = INFINITY;
		for (int i = 1; i < numObjs; i++)
		{
			newBBox     = objs[i]->getAABB();
			tempBBox    = AABB(tempBBox, newBBox);
			leftArea[i] = tempBBox.getArea();
		}

		// sweep from right
		tempBBox = AABB();
		rightArea[numObjs-1] = INFINITY;
		for (int i = numObjs-2; i >= 0; i--)
		{
			newBBox      = objs[i]->getAABB();
			tempBBox     = AABB(tempBBox, newBBox);
			rightArea[i] = tempBBox.getArea();
			thisCost     = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);
			if (thisCost < bestCost)
			{
				bestCost = thisCost;
				partPt   = i;
				bestAxis = 2;
			}
		}

		switch (bestAxis)
		{
		case 0 :
			qsort(objs, numObjs, sizeof(Object*), Object::sortByXComponent);
			break;
		case 1 :
			qsort(objs, numObjs, sizeof(Object*), Object::sortByYComponent);
			break;
		}
	}
#endif
}

void BVH_Node::buildSAH(Object** objs, u_int numObjs, float* leftArea, float* rightArea)
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
		nodeMin.x = min(nodeMin.x, objs[i]->getAABB().bbMin.x);
		nodeMin.y = min(nodeMin.y, objs[i]->getAABB().bbMin.y);
		nodeMin.z = min(nodeMin.z, objs[i]->getAABB().bbMin.z);
		nodeMax.x = max(nodeMax.x, objs[i]->getAABB().bbMax.x);
		nodeMax.y = max(nodeMax.y, objs[i]->getAABB().bbMax.y);
		nodeMax.z = max(nodeMax.z, objs[i]->getAABB().bbMax.z);
	}
	*bBox = AABB(nodeMin, nodeMax);

	u_int partPt = 0;
	if (numObjs <= MAX_LEAF_SIZE)		// Make a leaf...
	{
		numChildren = PUT_NUMCHILD(numObjs);
		isLeaf |= true;
		nodeCount++;
		leafCount++;

		this->objs = objs;
		qsort(objs, numObjs, sizeof(Object*), Object::sortByArea);
	}
	else
	{
		static int curDepth = 0;

		curDepth++;
		if (curDepth > maxDepth) maxDepth = curDepth; 
		nodeCount += 2;
		partitionSweepSAH(objs, numObjs, partPt, leftArea, rightArea);		// Find the best place to split the node

		numChildren = PUT_NUMCHILD(NODE_SIZE);										// Check out BVH.h, the first of numChildren is used by isLeaf.
		isLeaf |= false;

		Children = (BVH_Node*)_aligned_malloc(sizeof(BVH_Node)*2, 16);

		Object** leftObjs = objs;				// Get a new pointer so we don't have to do so much pointer arithmetic inside partitionSweep
		u_int leftNum = partPt+1;				// [0 .. partPt]

		Object** rightObjs = objs+(partPt+1);	// Same as before
		u_int rightNum = numObjs-partPt-1;		// [partPt+1 .. numObjs-1]

		Children[0].buildSAH(leftObjs, leftNum, leftArea, rightArea);	// Build left and right child nodes.
		Children[1].buildSAH(rightObjs, rightNum, leftArea, rightArea);
		curDepth--;
		return;
	}
}

void BVH_Node::partitionSweepSAH(Object** objs, u_int numObjs, u_int& partPt, float* leftArea, float* rightArea)
{
	float bestCost  = MIRO_TMAX;
	float thisCost  = 0;
	int bestAxis	= -1;

	// Using ideas from Wald's "Ray Tracing Deformable Scenes Using Dynamic Bounding Volume Hierarchies" (Siggraph, 2007).

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
		thisCost = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);		// Calculate the cost of this partition.
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
		thisCost = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);
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
		thisCost = calcSAHCost( (i+1), leftArea[i], (numObjs-i-1), rightArea[i]);
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

float BVH_Node::calcSAHCost(int leftNum, float leftArea, int rightNum, float rightArea)
{
#ifdef USE_TRI_PACKETS	// Try to force full 4-leaves if using SIMD to intersect 4 tris at same time.
						// These values have to be tested, it's very dependent on scene structure. 
						// Current values are pretty good for Sponza using full SAH build (works OK for binned build also)
	if (leftNum + rightNum >= 256)
	{
		return ((float)leftNum)*leftArea + ((float)rightNum)*rightArea;
	}

	else if (leftNum + rightNum >= 64)
	{
		float penaltyMult = 1.0;
		if (leftNum % 16 != 0 && rightNum % 16 != 0) penaltyMult = 1.25;	// Try to maintain multiple of 16 nodes by (slightly) penalizing other node sizes
		return penaltyMult*(((float)leftNum)*leftArea + ((float)rightNum)*rightArea);
	}

	float penaltyMult = 1.0;
	if (leftNum % 16 != 0 && rightNum % 16 != 0) penaltyMult = 1.25;	// Try to maintain multiple of 16 nodes by (slightly) penalizing other node sizes
	if (leftNum % 4 != 0 && rightNum % 4 != 0) penaltyMult *= 2;		// Try to maintain multiple of 4 leaves by penalizing other node sizes (slightly more)

	return penaltyMult*(((float)leftNum)*leftArea + ((float)rightNum)*rightArea);

#else

	return ((float)leftNum)*leftArea + ((float)rightNum)*rightArea;

#endif
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

#ifdef USE_TRI_PACKETS
bool intersect4(HitInfo& result, const Ray& r, float tMin, BVH_Node::TriCache4* triCache);
#endif

bool BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin) const 
{
	if (!use_BVH)
	{
		bool hit = false;

		for (size_t i = 0; i < m_objects->size(); ++i)
		{
			if ((*m_objects)[i]->intersect(minHit, ray, tMin))
			{
				hit = true;
			}
		}    
		return hit;
	}
	else
	{
		bool hit = false;
		int stackIndex = -1;
	    BVH_Node* BVH_Stack[256];

		// Push children on "stack"
		for (int i = NODE_SIZE-1; i >= 0; i--)
		{		
			++stackIndex;
			BVH_Stack[stackIndex] = &m_baseNode.Children[i];
		}

		bool isFirst = true;
		while (stackIndex >= 0)
		{
			// Intersect first child box

			//#pragma omp atomic
			BVH::rayBoxIntersections++;

			if (BVH_Stack[stackIndex]->intersect(minHit, ray, tMin))
			{
				if (BVH_Stack[stackIndex]->isLeaf & true)		// If this is a leaf, then we need to intersect the objects inside...
				{
#ifdef USE_TRI_PACKETS
					if (intersect4(minHit, ray, tMin, BVH_Stack[stackIndex]->triCache))
					{
						hit = true;
					}
#else
					for (int i = 0; i < GET_NUMCHILD(BVH_Stack[stackIndex]->numChildren); i++)
					{
						if (BVH_Stack[stackIndex]->objs[i]->intersect(minHit, ray, tMin))
						{
							hit = true;
						}
					}
#endif
					stackIndex--;
				}
				else
				{
					// Push on stack by overwriting "this" node (same as popping this node before pushing)
					// I go from last to first - otherwise I'd overwrite the current pointer first and results are wrong.
					for (int i = NODE_SIZE-1; i >= 0; i--)
					{				
						BVH_Stack[stackIndex+i] = &BVH_Stack[stackIndex]->Children[NODE_SIZE-1-i];
					}

					stackIndex += NODE_SIZE-1;
				}
			}
			else 
			{
				stackIndex--;
			}
		}
		return hit;
	}
}

bool BVH_Node::intersect(HitInfo& result, const Ray& ray, float tMin) const
{
/*#ifndef NO_SSE												// Implementation of Manny Ko's algorithm, from http://tog.acm.org/resources/RTNews/html/rtnv23n1.html#art7
																// Had to guess as to a few details, but works great.
																// THIS DOESN'T WORK ANYMORE AND CAN'T FIGURE OUT WHY - it actually works but slows everything down to a crawl
	const __m128 minPoint  = loadps(&bBox->bbMin.x);
	const __m128 maxPoint  = loadps(&bBox->bbMax.x);
	
	const __m128 t0v = mulps(subps(minPoint, ray._o), ray._id);
	const __m128 t1v = mulps(subps(maxPoint, ray._o), ray._id);

	const __m128 interval_min = maxss(hMax(minps(t0v, t1v)), loadss(&tMin)); // tMin
	const __m128 interval_max = minss(hMin(maxps(t0v, t1v)), loadss(&result.t)); // tMax

	return (_mm_comile_ss(interval_min, interval_max));

#else*/															// Standard slabs test, following Williams et al. http://portal.acm.org/citation.cfm?id=1198748
	
	float bbmin, bbmax, tymin, tymax, tzmin, tzmax;	

	bbmin = ((&bBox->bbMin)[ray.bounces_flags & SIGN_X_FLAG].x - ray.o[0]) * ray.id[0];
	bbmax = ((&bBox->bbMin)[1-(ray.bounces_flags & SIGN_X_FLAG)].x - ray.o[0]) * ray.id[0];
	tymin = ((&bBox->bbMin)[GET_Y_SIGN(ray.bounces_flags & SIGN_Y_FLAG)].y - ray.o[1]) * ray.id[1];
	tymax = ((&bBox->bbMin)[1-(GET_Y_SIGN(ray.bounces_flags & SIGN_Y_FLAG))].y - ray.o[1]) * ray.id[1];

	if ((bbmin > tymax) || (tymin > bbmax))
		return false;

	if (tymin > bbmin)
		bbmin = tymin;
	if (tymax < bbmax)
		bbmax = tymax;

	tzmin = ((&bBox->bbMin)[GET_Z_SIGN(ray.bounces_flags & SIGN_Z_FLAG)].z - ray.o[2]) * ray.id[2];
	tzmax = ((&bBox->bbMin)[1-(GET_Z_SIGN(ray.bounces_flags & SIGN_Z_FLAG))].z - ray.o[2]) * ray.id[2];

	if ((bbmin > tzmax) || (tzmin > bbmax))
		return false;
	if (tzmin > bbmin)
		bbmin = tzmin;
	if (tzmax < bbmax)
		bbmax = tzmax;

	return ((bbmin <= result.t) && (bbmax >= tMin));
//#endif
}

#ifdef USE_TRI_PACKETS
bool intersect4(HitInfo& result, const Ray& r, float tMin, BVH_Node::TriCache4* triCache)
{
//#pragma omp atomic
	Ray::rayTriangleIntersections++;

	bool shadow = r.bounces_flags & IS_SHADOW_RAY;

	// Cross product
	const __m128 pVecx4		= subps(mulps(r.dy4, triCache->edge1z4), mulps(r.dz4, triCache->edge1y4)); 
	const __m128 pVecy4		= mulps(_negonesps, subps(mulps(r.dx4, triCache->edge1z4), mulps(r.dz4, triCache->edge1x4)));
	const __m128 pVecz4		= subps(mulps(r.dx4, triCache->edge1y4), mulps(r.dy4, triCache->edge1x4));

	const __m128 det4		= SoADot(triCache->edge0x4, triCache->edge0y4, triCache->edge0z4, pVecx4, pVecy4, pVecz4);
	const __m128 invDet4	= recipps(det4);

	const __m128 tVecx4		= subps(r.ox4, triCache->Ax4);
	const __m128 tVecy4		= subps(r.oy4, triCache->Ay4);
	const __m128 tVecz4		= subps(r.oz4, triCache->Az4);

	const __m128 a4			= mulps(invDet4, SoADot(tVecx4, tVecy4, tVecz4, pVecx4, pVecy4, pVecz4));
	int aMask;
	if ((aMask = (movemaskps(cmpgteqps(a4, _zerosps)) & movemaskps(cmpleqps(a4, _onesps)))) == 0) return false;

	// Cross product
	const __m128 qVecx4		= subps(mulps(tVecy4, triCache->edge0z4), mulps(tVecz4, triCache->edge0y4)); 
	const __m128 qVecy4		= mulps(_negonesps, subps(mulps(tVecx4, triCache->edge0z4), mulps(tVecz4, triCache->edge0x4)));
	const __m128 qVecz4		= subps(mulps(tVecx4, triCache->edge0y4), mulps(tVecy4, triCache->edge0x4));

	const __m128 b4 = mulps(invDet4, SoADot(r.dx4, r.dy4, r.dz4, qVecx4, qVecy4, qVecz4));
	int bMask;
	if ((bMask = (aMask & movemaskps(cmpgteqps(b4, _zerosps)) & movemaskps(cmpleqps(b4, _onesps)) & movemaskps(cmpleqps(addps(a4, b4), _onesps)))) == 0) return false;

	const __m128 tMin4  = setSSE(tMin);
	const __m128 tMax4  = setSSE(result.t);

	const __m128 newT4 = mulps(invDet4, SoADot(triCache->edge1x4, triCache->edge1y4, triCache->edge1z4, qVecx4, qVecy4, qVecz4));
	int tMask;
	if ((tMask = (bMask & movemaskps(cmpgteqps(newT4, tMin4)) & movemaskps(cmplessps(newT4, tMax4)))) == 0) return false;

	if (shadow) return true;

	ALIGN_SSE float newT[4];
	ALIGN_SSE float a[4];
	ALIGN_SSE float b[4];
	storeps(newT4, newT);
	storeps(a4, a);
	storeps(b4, b);

	for (int i = 0; i < 4; i++)
	{
		newT[i] = (tMask & 1<<i) ? newT[i] : MIRO_TMAX;
	}

	float lowest = newT[0];
	int lindex = 0;
	for (int i = 1; i < 4; i++ ) if (newT[i] < lowest)
	{
		lowest = newT[i];
		lindex = i;
	}

	result.t	=			newT[lindex];
	result.a	=			   a[lindex];
	result.b	=			   b[lindex];
	result.obj	= triCache->tris[lindex];

	return true;
}
#endif