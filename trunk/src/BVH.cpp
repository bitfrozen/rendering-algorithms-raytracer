#include "BVH.h"
#include "Ray.h"
#include "Console.h"
#include "TriangleMesh.h"
#include "Triangle.h"

void
BVH::build(Objects * objs)
{
    m_objects = objs;

	Vector3 sceneMin = Vector3(MIRO_TMAX);
	Vector3 sceneMax = Vector3(-MIRO_TMAX);	
	Vector3 tempMin, tempMax;

	Objects::const_iterator objIter;
	for (objIter = objs->begin(); objIter != objs->end(); objIter++)
	{
		(*objIter)->getBounds(tempMin, tempMax);
		sceneMin.x = std::min(sceneMin.x, tempMin.x);
		sceneMin.y = std::min(sceneMin.y, tempMin.y);
		sceneMin.z = std::min(sceneMin.z, tempMin.z);
		sceneMax.x = std::max(sceneMax.x, tempMax.x);
		sceneMax.y = std::max(sceneMax.y, tempMax.y);
		sceneMax.z = std::max(sceneMax.z, tempMax.z);
	}
	m_baseNode.vMin = sceneMin;
	m_baseNode.vMax = sceneMax;
	m_baseNode.isLeaf = false;
	m_baseNode.obj = NULL;

	for (objIter = objs->begin(); objIter != objs->end(); objIter++)
	{
		BVH_Node* newChild = new BVH_Node();
		m_baseNode.Children.push_back(*newChild);
		m_baseNode.Children.back().build(*objIter);
	}	
}

bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    /*// Here you would need to traverse the BVH to perform ray-intersection
    // acceleration. For now we just intersect every object.
    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = tempMinHit.t = tMax;
    
    for (size_t i = 0; i < m_objects->size(); ++i)
    {
        if ((*m_objects)[i]->intersect(tempMinHit, ray, tMin, tMax))
        {
            hit = true;
            if (tempMinHit.t < minHit.t)
                minHit = tempMinHit;
        }
    }    
    return hit;*/
	return m_baseNode.intersect(minHit, ray, tMin, tMax);
}

void BVH_Node::build(Object* obj, int first /* = 0 */, int last /* = -100 */)
{
	Triangle* triMesh;
	if (triMesh = dynamic_cast<Triangle*>(obj))
	{
		if (last == -100)  // First run, build AABB for whole object
		{
			triMesh->getBounds(vMin, vMax);
			int numFaces = triMesh->getNumTris();
			isLeaf = false;

			BVH_Node* newChild = new BVH_Node();
			Children.push_back(*newChild);
			Children.back().build(triMesh, 0, numFaces/2-1);

			newChild = new BVH_Node();
			Children.push_back(*newChild);
			Children.back().build(triMesh, numFaces/2, numFaces-1);

			this->obj = NULL;
		} 
		else if (last-first > 0)
		{
			triMesh->getBounds(vMin, vMax, first, last);
			int numFaces = last-first+1;
			isLeaf = false;

			BVH_Node* newChild = new BVH_Node();
			Children.push_back(*newChild);
			Children.back().build(triMesh, first, first+numFaces/2-1);

			newChild = new BVH_Node();
			Children.push_back(*newChild);
			Children.back().build(triMesh, first+numFaces/2, last);

			this->obj = NULL;
		}
		else if (first == last)
		{
			triMesh->getBounds(vMin, vMax, first, last);
			isLeaf = true;
			this->obj = new Triangle(triMesh->getMesh(), first);
			this->obj->setMaterial(triMesh->getMaterial());
		}		
	} 
	else
	{
		obj->getBounds(vMin, vMax);
		isLeaf = true;
		this->obj = obj;
	}
	
}

bool BVH_Node::intersect(HitInfo& result, const Ray& ray, float tMin, float tMax) const
{
	float bbmin, bbmax, tymin, tymax, tzmin, tzmax;
	if (ray.d.x >= 0)
	{
		bbmin = (vMin.x - ray.o.x) / ray.d.x;
		bbmax = (vMax.x - ray.o.x) / ray.d.x;
	}
	else
	{
		bbmin = (vMax.x - ray.o.x) / ray.d.x;
		bbmax = (vMin.x - ray.o.x) / ray.d.x;
	}
	if (ray.d.y >= 0)
	{
		tymin = (vMin.y - ray.o.y) / ray.d.y;
		tymax = (vMax.y - ray.o.y) / ray.d.y;
	}
	else
	{
		tymin = (vMax.y - ray.o.y) / ray.d.y;
		tymax = (vMin.y - ray.o.y) / ray.d.y;
	}
	if ((bbmin > tymax) || (tymin > bbmax))
		return false;
	
	if (tymin > bbmin)
		bbmin = tymin;
	if (tymax < bbmax)
		bbmax = tymax;

	if (ray.d.z >= 0)
	{
		tzmin = (vMin.z - ray.o.z) / ray.d.z;
		tzmax = (vMax.z - ray.o.z) / ray.d.z;
	}
	else
	{
		tzmin = (vMax.z - ray.o.z) / ray.d.z;
		tzmax = (vMin.z - ray.o.z) / ray.d.z;
	}
	if ((bbmin > tzmax) || (tzmin > bbmax))
		return false;

	if (tzmin > bbmin)
		bbmin = tzmin;
	if (tzmax < bbmax)
		bbmax = tzmax;

	bool hit = ((bbmin < tMax) && (bbmax > tMin));

	if (hit)
	{
		if (isLeaf)
		{
			return obj->intersect(result, ray, tMin, tMax);
		}
		else
		{
			hit = false;
			HitInfo tempMinHit;
			tempMinHit.t = tMax;

			std::vector<BVH_Node>::const_iterator childIter;
			for (childIter = Children.begin(); childIter !=Children.end(); childIter++)
			{				
				if (childIter->intersect(tempMinHit, ray, tMin, tMax))
				{
					hit = true;
					if (tempMinHit.t < result.t)
						result = tempMinHit;
				}				
			}
			return hit;
		}
	}
	return false;
}