#ifndef CSE168_RAY_H_INCLUDED
#define CSE168_RAY_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Material.h"

// Bitfield flags
#define SIGN_X_FLAG			0x1
#define GET_X_SIGN(a)		(a)>>0
#define SIGN_Y_FLAG			0x2
#define GET_Y_SIGN(a)		(a)>>1
#define SIGN_Z_FLAG			0x4
#define GET_Z_SIGN(a)		(a)>>2
#define IS_REFLECT_RAY		0x8
#define IS_REFRACT_RAY		0x10
#define IS_PRIMARY_RAY		0x20
#define IS_GI_RAY			0x40
#define BOUNCES_MASK		0x7F80		// 8-bit mask, up to 255 bounces
#define GET_BOUNCES(a)		(a)>>7
#define IS_SHADOW_RAY		0x8000

ALIGN_SSE class Ray
{
public:
    float ox, oy, oz, ow;      //!< Origin of ray
	float dx, dy, dz, dw;      //!< Direction of ray
	float idx, idy, idz, idw;  // Reciprocal of direction, used for slabs test (AABB intersection)
	unsigned int bounces_flags;
	typedef std::vector<float> IORList;		
	IORList r_IOR;							// History of IOR this ray has traversed..

    Ray()
    {
		ox = oy = oz = 0.0f; ow = 1.0f;
		dx = dy = dw = 0.0f; dz = 1.0f;
		idx = idy = idw = 0.0f; idz = 1.0f;
		bounces_flags = IS_PRIMARY_RAY;
		r_IOR.push_back(1.001f);
    }

    Ray(const Vector3& o, const Vector3& d, float IOR = 1.001f, unsigned int bounces = 0, unsigned int flags = IS_PRIMARY_RAY)
    {
		ox = o.x; oy = o.y; oz = o.z; ow = 1.0f;
		dx = d.x; dy = d.y; dz = d.z; dw = 0.0f;
		idx = 1.0f/dx; idy = 1.0f/dy; idz = 1.0f/dz; idw = 0.0f;
		bounces_flags = bounces<<7 | (idz < 0)<<2 | (idy < 0)<< 1 | (idx < 0) | flags;
		if (!(flags & IS_SHADOW_RAY))
			r_IOR.push_back(IOR);
    }

	Ray(const Vector3& o, const Vector3& d, const IORList& IOR, unsigned int bounces = 0, unsigned int flags = IS_PRIMARY_RAY)
	{
		ox = o.x; oy = o.y; oz = o.z; ow = 1.0f;
		dx = d.x; dy = d.y; dz = d.z; dw = 0.0f;
		idx = 1.0f/dx; idy = 1.0f/dy; idz = 1.0f/dz; idw = 0.0f;
		bounces_flags = bounces<<7 | (idz < 0)<<2 | (idy < 0)<< 1 | (idx < 0) | flags;
		r_IOR = IOR;
	}
};


//! Contains information about a ray hit with a surface.
/*!
    HitInfos are used by object intersection routines. They are useful in
    order to return more than just the hit distance.
*/
ALIGN_SSE class HitInfo
{
public:
	float px, py, pz, pw;
    float t, a, b;                      //!< The hit distance, u and v barycentric coordinates
    Vector3 N, geoN;                    //!< Shading(interpolated) and geometric normal vector
    const Material* material;           //!< Material of the intersected object

    //! Default constructor.
    explicit HitInfo(float t = MIRO_TMAX, float a = 0.0f, float b = 0.0f,
                     float px = 0.0f, float py = 0.0f, float pz = 0.0f, float pw = 0.0f,
                     const Vector3& N = Vector3(0.0f, 1.0f, 0.0f)) :
        t(t), a(a), b(b), px(px), py(py), pz(pz), pw(pw), N(N), geoN(N), material (0)
    {
        // empty
    }
};

#endif // CSE168_RAY_H_INCLUDED
