#include "Ray.h"
#include "ProxyObject.h"
#include "TriangleMesh.h"

const void HitInfo::getAllInfos(Vector3 &N, Vector3 &geoN, Vector3 &T, Vector3 &BT, float &uCoord, float &vCoord) const
{
	if (!obj)												// Do nothing if there's no object
	{
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;

	ti3           = theMesh->m_vertexIndices[meshIndex];	// Get the geometric normal
	Vector3 edge0 = theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x];
	Vector3 edge1 = theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x];
	geoN		  = cross(edge0, edge1).normalized();
	geoN.__dummy  = 0.f;

	ti3       = theMesh->m_normalIndices[meshIndex];		// Get the interpolated normal
	float c   = 1.0f-a-b;
	N		  = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b).normalized());
	N.__dummy = 0.f;

	if (m_proxy)
	{
		geoN = (m_proxy->getMatrix().m_invTranspose * geoN).normalized();
		N = (m_proxy->getMatrix().m_invTranspose * N).normalized();
	}

	if (theMesh->m_texCoordIndices)							// If possible, get the interpolated u, v coordinates
	{
		T		  = Vector3((theMesh->m_tangents[ti3.x]*c+theMesh->m_tangents[ti3.y]*a+theMesh->m_tangents[ti3.z]*b).normalized());
		BT		  = Vector3((theMesh->m_biTangents[ti3.x]*c+theMesh->m_biTangents[ti3.y]*a+theMesh->m_biTangents[ti3.z]*b).normalized());

		ti3    = theMesh->m_texCoordIndices[meshIndex];
		uCoord = theMesh->m_texCoords[ti3.x].x*c+theMesh->m_texCoords[ti3.y].x*a+theMesh->m_texCoords[ti3.z].x*b;
		vCoord = theMesh->m_texCoords[ti3.x].y*c+theMesh->m_texCoords[ti3.y].y*a+theMesh->m_texCoords[ti3.z].y*b;
	}
	else													// We always return texture coordinates
	{
		T = 0;
		BT = 0;
		uCoord = a;
		vCoord = b;
	}
}

const void HitInfo::getInterpolatedNormal(Vector3& N) const
{
	if (!obj)												// Do nothing if there's no object
	{
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;

	ti3       = theMesh->m_normalIndices[meshIndex];		// Get the interpolated normal
	float c   = 1.0f-a-b;
	N		  = Vector3((theMesh->m_normals[ti3.x]*c+theMesh->m_normals[ti3.y]*a+theMesh->m_normals[ti3.z]*b).normalized());
}

const void HitInfo::getInterpolatedTangent(Vector3& T) const
{
	if (!obj || !obj->m_mesh->m_texCoordIndices)			// Do nothing if there's no object or no texture coords
	{
		T = 0;
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;

	ti3       = theMesh->m_normalIndices[meshIndex];		// Get the interpolated normal
	float c   = 1.0f-a-b;
	T		  = Vector3((theMesh->m_tangents[ti3.x]*c+theMesh->m_tangents[ti3.y]*a+theMesh->m_tangents[ti3.z]*b).normalized());
}

const void HitInfo::getInterpolatedBiTangent(Vector3& BT) const
{
	if (!obj || !obj->m_mesh->m_texCoordIndices)			// Do nothing if there's no object or no texture coords
	{
		BT = 0;
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;

	ti3       = theMesh->m_normalIndices[meshIndex];		// Get the interpolated normal
	float c   = 1.0f-a-b;
	BT		  = Vector3((theMesh->m_biTangents[ti3.x]*c+theMesh->m_biTangents[ti3.y]*a+theMesh->m_biTangents[ti3.z]*b).normalized());
}

const void HitInfo::getGeoNormal(Vector3& geoN) const
{
	if (!obj)												// Do nothing if there's no object
	{
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;

	ti3           = theMesh->m_vertexIndices[meshIndex];	// Get the geometric normal
	Vector3 edge0 = theMesh->m_vertices[ti3.y] - theMesh->m_vertices[ti3.x];
	Vector3 edge1 = theMesh->m_vertices[ti3.z] - theMesh->m_vertices[ti3.x];
	geoN		  = cross(edge0, edge1).normalized();
}

const void HitInfo::getUVs(float& uCoord, float &vCoord) const
{
	if (!obj)												// Do nothing if there's no object
	{
		return;
	}

	TriangleMesh::TupleI3 ti3;								// Get pointer and index for the mesh
	TriangleMesh* theMesh = obj->m_mesh;
	u_int meshIndex       = obj->m_index;
	float c   = 1.0f-a-b;

	if (theMesh->m_texCoordIndices)							// If possible, get the interpolated u, v coordinates
	{
		ti3    = theMesh->m_texCoordIndices[meshIndex];
		uCoord = theMesh->m_texCoords[ti3.x].x*c+theMesh->m_texCoords[ti3.y].x*a+theMesh->m_texCoords[ti3.z].x*b;
		vCoord = theMesh->m_texCoords[ti3.x].y*c+theMesh->m_texCoords[ti3.y].y*a+theMesh->m_texCoords[ti3.z].y*b;
	}
	else													// We always return texture coordinates
	{
		uCoord = a;
		vCoord = b;
	}
}