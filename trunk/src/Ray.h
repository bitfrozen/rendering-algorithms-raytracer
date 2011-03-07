#ifndef CSE168_RAY_H_INCLUDED
#define CSE168_RAY_H_INCLUDED

#include "ProxyObject.h"
#include <vector>
#include "Vector3.h"
#include "SSE.h"
#include <omp.h>

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
#define GET_BOUNCES(a)		(a & BOUNCES_MASK)>>7
#define IS_SHADOW_RAY		0x8000
#define GI_BOUNCES_MASK		0xFF0000	// 8-bit mask, up to 255 bounces
#define GET_GI_BOUNCES(a)	(a & GI_BOUNCES_MASK)>>16

ALIGN_64 class Ray
{
public:
	union {float  o[4]; __m128  _o;};	//!< Origin of ray
	union {float  d[4]; __m128  _d;};	//!< Direction of ray
	union {float id[4]; __m128 _id;};	// Reciprocal of direction, used for slabs test (AABB intersection)
#ifdef USE_TRI_PACKETS
	__m128  ox4,  oy4,  oz4;
	__m128  dx4,  dy4,  dz4;
	__m128 idx4, idy4, idz4;
#endif
	static unsigned long int counter[256], rayTriangleIntersections[256];
	unsigned int bounces_flags;
	float time;
	
	struct IORList {
		float _IORList[256];
		unsigned int IOR_idx;
		IORList() {_IORList[0] = 1.0f; IOR_idx = 0;}
		float operator()() {return _IORList[IOR_idx];}
		void pop() {if (IOR_idx > 0) IOR_idx--;}
		void push(float num) {_IORList[++IOR_idx] = num;}
	};
	mutable IORList r_IOR;							// History of IOR this ray has traversed..

    Ray(const unsigned int threadID = 0, const float t = 0.f)
    {
		//#pragma omp atomic
		counter[threadID]++;

		 o[0] =  o[1] =  o[2] = o[3] = 0.0f;
		 d[0] =  d[1] =  d[3] = 0.0f;  d[2] = 1.0f;
		id[0] = id[1] = id[3] = 0.0f; id[2] = 1.0f;
#ifdef USE_TRI_PACKETS
		 ox4 = setZero;  oy4 = setZero;  oz4 = setZero;
		 dx4 = setZero;  dy4 = setZero;  dz4 = setZero;
		idx4 = setZero; idy4 = setZero; idz4 = setZero;
#endif
		bounces_flags = IS_PRIMARY_RAY;
		r_IOR.push(1.001f);
		time = t;
    }

    Ray(const unsigned int threadID, const Vector3& o, const Vector3& d, const float t = 0.f, float IOR = 1.001f, unsigned int bounces = 0, unsigned int giBounces = 0, unsigned int flags = IS_PRIMARY_RAY)
    {
		//#pragma omp atomic
		counter[threadID]++;

		this->o[0] = o.x; this->o[1] = o.y; this->o[2] = o.z; this->o[3] = 1.0f;
		this->d[0] = d.x; this->d[1] = d.y; this->d[2] = d.z; this->d[3] = 0.0f;

		id[0] = 1.0f/d[0]; if (d[0] == 0.f)
		{
			id[0] = (id[0] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[1] = 1.0f/d[1]; if (d[1] == 0.f)
		{
			id[1] = (id[1] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[2] = 1.0f/d[2]; if (d[2] == 0.f)
		{
			id[2] = (id[2] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[3] = 1.0f;
#ifdef USE_TRI_PACKETS
		 ox4 = setSSE( o[0]);  oy4 = setSSE( o[1]);  oz4 = setSSE( o[2]);
		 dx4 = setSSE( d[0]);  dy4 = setSSE (d[1]);  dz4 = setSSE( d[2]);
		idx4 = setSSE(id[0]); idy4 = setSSE(id[1]); idz4 = setSSE(id[2]);
#endif
		bounces_flags = giBounces<<16 | bounces<<7 | (id[2] < 0)<<2 | (id[1] < 0)<< 1 | (id[0] < 0) | flags;
		if (!(flags & IS_SHADOW_RAY))
			r_IOR.push(IOR);
		time = t;
    }

	Ray(const unsigned int threadID, const Vector3& o, const Vector3& d, const float t, const IORList IOR, unsigned int bounces = 0, unsigned int giBounces = 0, unsigned int flags = IS_PRIMARY_RAY)
	{
		//#pragma omp atomic
		r_IOR = IOR;
		counter[threadID]++;

		this->o[0] = o.x; this->o[1] = o.y; this->o[2] = o.z; this->o[3] = 1.0f;
		this->d[0] = d.x; this->d[1] = d.y; this->d[2] = d.z; this->d[3] = 0.0f;

		id[0] = 1.0f/d[0]; if (d[0] == 0.f)
		{
			id[0] = (id[0] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[1] = 1.0f/d[1]; if (d[1] == 0.f)
		{
			id[1] = (id[1] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[2] = 1.0f/d[2]; if (d[2] == 0.f)
		{
			id[2] = (id[2] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[3] = 1.0f;

#ifdef USE_TRI_PACKETS
		 ox4 = setSSE( o[0]);  oy4 = setSSE( o[1]);  oz4 = setSSE( o[2]);
		 dx4 = setSSE( d[0]);  dy4 = setSSE (d[1]);  dz4 = setSSE( d[2]);
		idx4 = setSSE(id[0]); idy4 = setSSE(id[1]); idz4 = setSSE(id[2]);
#endif
		bounces_flags = giBounces<<16 | bounces<<7 | (id[2] < 0)<<2 | (id[1] < 0)<< 1 | (id[0] < 0) | flags;
		time = t;
	}

	__forceinline void set(const unsigned int threadID, const Vector3& o, const Vector3& d, const float t = 0.f, const float IOR = 1.001f, unsigned int bounces = 0, unsigned int giBounces = 0, unsigned int flags = IS_PRIMARY_RAY)
	{
		//#pragma omp atomic
		counter[threadID]++;

		this->o[0] = o.x; this->o[1] = o.y; this->o[2] = o.z; this->o[3] = 1.0f;
		this->d[0] = d.x; this->d[1] = d.y; this->d[2] = d.z; this->d[3] = 0.0f;

		id[0] = 1.0f/d[0]; if (d[0] == 0.f)
		{
			id[0] = (id[0] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[1] = 1.0f/d[1]; if (d[1] == 0.f)
		{
			id[1] = (id[1] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[2] = 1.0f/d[2]; if (d[2] == 0.f)
		{
			id[2] = (id[2] < -0.f) ? -MIRO_TMAX : MIRO_TMAX;
		}
		id[3] = 1.0f;

#ifdef USE_TRI_PACKETS
		ox4 = setSSE( o[0]);  oy4 = setSSE( o[1]);  oz4 = setSSE( o[2]);
		dx4 = setSSE( d[0]);  dy4 = setSSE (d[1]);  dz4 = setSSE( d[2]);
		idx4 = setSSE(id[0]); idy4 = setSSE(id[1]); idz4 = setSSE(id[2]);
#endif
		bounces_flags = giBounces<<16 | bounces<<7 | (id[2] < 0)<<2 | (id[1] < 0)<< 1 | (id[0] < 0) | flags;
		r_IOR.pop();
		r_IOR.push(IOR);
		time = t;
	}

	__forceinline const Vector3 getPoint(const float t) const
	{
		Vector3 P;
#ifndef NO_SSE
		storeps(addps(_o, mulps(setSSE(t), _d)), P.v);	// Hit position
#else
		P = Vector3(o[0] + t*d[0], o[1] + t*d[1], o[2] + t*d[2]);
#endif
		return P;
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
	Object* obj;						//!< Pointer to intersected object
	ProxyObject* m_proxy;
    float t, a, b;                      //!< The hit distance, a and b barycentric coordinates

    //! Default constructor.
    explicit HitInfo(float t = MIRO_TMAX, float a = 0.0f, float b = 0.0f, Object* obj = NULL, ProxyObject* proxy = NULL) : t(t), a(a), b(b), obj(obj), m_proxy(proxy){};
	const void getAllInfos(Vector3 &N, Vector3 &geoN, Vector3 &T, Vector3 &BT, float &uCoord, float &vCoord) const;
	const void getInterpolatedNormal(Vector3& N) const;
	const void getInterpolatedTangent(Vector3& T) const;
	const void getInterpolatedBiTangent(Vector3& BT) const;
	const void getGeoNormal(Vector3& geoN) const;
	const void getUVs(float& uCoord, float &vCood) const;
};

#endif // CSE168_RAY_H_INCLUDED