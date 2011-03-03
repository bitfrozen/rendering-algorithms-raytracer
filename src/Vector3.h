#ifndef CSE168_VECTOR3_H_INCLUDED
#define CSE168_VECTOR3_H_INCLUDED

#include <math.h>
#include <float.h>
#include <iostream>
#include "Miro.h"
#include "SSE.h"

#ifdef WIN32
#pragma warning(disable:4305) // disable useless warnings
#pragma warning(disable:4244)
#endif

class Vector3
{

public:
	union {float v[4]; __m128 _v; struct {float x, y, z, __dummy;}; };      // The x & y & z coordinates.

    Vector3() :
        x(0), y(0), z(0), __dummy(1) {}

    Vector3(float s) :
        x(s), y(s), z(s), __dummy(1) {}

    Vector3(float xVal, float yVal, float zVal) :
        x(xVal), y(yVal), z(zVal), __dummy(1) {}

	Vector3(const __m128& vec) : _v(vec) {}

    //! Assignment operator.
    /*!
        Assigns the values from \a a to this Vec3.
    */
    __forceinline const Vector3 & operator=(const Vector3& a) {x = a.x; y = a.y; z = a.z; return *this;}
    
    //! Assignment operator.
    /*!
        Sets all components of this Vec3 to \a a.
    */
    __forceinline const Vector3 & operator=(float a) {x = y = z = a; return *this;}

    __forceinline void set(float a) {x = y = z = a;}
    __forceinline void set(float a, float b, float c, float d = 0.0) {x = a; y = b; z = c;}
    __forceinline void set(const Vector3 & a) {x = a.x; y = a.y; z = a.z;}
    
    
    //! Access operator.        
    /*!
        Returns the ith component of the vector.
        \param i The component to return.
        \warning i must be either 0, 1, or 2 in order to get expected results.
    */
    __forceinline float & operator[](int i) {return v[i];}
    
    //! Constant access operator.
    /*!
        Returns the ith component of a constant vector.
        \param i The component to return.
        \warning i must be either 0, 1, or 2 in order to get expected results.
    */
    __forceinline const float & operator[](int i) const {return v[i];}


    //! Component-wise vector addition operator.
    __forceinline Vector3 operator+(const Vector3& a) const
    {
        return Vector3(x + a.x, y + a.y, z + a.z);
    }
    
    //! Component-wise vector addition-assignment operator.
    __forceinline const Vector3 & operator+=(const Vector3& a)
    {
        x += a.x; y += a.y; z += a.z; return *this;
    }

    //! Scalar addition-assignment operator.
    __forceinline const Vector3 & operator+=(float a) {x += a; y += a; z += a; return *this;}


    //! Component-wise vector subtraction operator.
    __forceinline Vector3 operator-(const Vector3& a) const
    {
        return Vector3(x - a.x, y - a.y, z - a.z);
    }
    
    //! Component-wise vector subtraction-assignment operator.
    __forceinline const Vector3 & operator-=(const Vector3& a)
    {
        x -= a.x; y -= a.y; z -= a.z; return *this;
    }
    
    //! Component-wise scalar subtraction assignment operator.
    __forceinline const Vector3 & operator-=(float a) {x -= a; y -= a; z -= a; return *this;}


    //! Scalar multiplication operator.
    __forceinline Vector3 operator*(float a) const {return Vector3(x * a, y * a, z * a);}
    
    //! Component-wise vector multiplication operator.
    __forceinline Vector3 operator*(const Vector3& a) const
    {
        return Vector3(x * a.x, y * a.y, z * a.z);
    }
    
    //! Scalar multiplication-assignment operator.
    __forceinline const Vector3 & operator*=(float a) {x *= a; y *= a; z *= a; return *this;}
    
    //! Component-wise vector multiplication-assignment operator.
    __forceinline const Vector3 & operator*=(const Vector3& a)
    {
        x *= a.x; y *= a.y; z *= a.z; return *this;
    }
    
    //! Negation operator.
    __forceinline Vector3 operator-() const {return Vector3(-x, -y, -z);}
    __forceinline const Vector3 & negate() {x = -x; y = -y; z = -z; return *this;}


    //! Scalar division operator.
    __forceinline Vector3 operator/(float a) const
    {
        float inv = float(1) / a;
        return Vector3(x * inv, y * inv, z * inv);
    }
    
    //! Component-wise vector division operator.
    __forceinline Vector3 operator/(const Vector3 & a) const
    {
        return Vector3(x / a.x, y / a.y, z / a.z);
    }
    
    //! Scalar division-assignment operator.
    __forceinline const Vector3 & operator/=(float a)
    {
        float inv = float(1) / a;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    
    //! Component-wise vector division-assignment operator.
    __forceinline const Vector3 & operator/=(const Vector3 & a)
    {
        x /= a.x; y /= a.y; z /= a.z; return *this;
    }


    //! Vector equivalence operator.
    /*!
        Tests to see if each component of \a v is equal to each component of
        this Vector3.
    */
    __forceinline bool operator==(const Vector3 & a) const
    {
        return(a.x == x && a.y == y && a.z == z);
    }
    
    //! Vector difference operator.
    /*!
        Tests to see if any component is different between the two Vec3s.
    */
    __forceinline bool operator!=(const Vector3 & a) const
    {
        return(a.x != x || a.y != y || a.z != z);
    }


    //! Length<sup>2</sup>.
    /*!
        Returns the geometric length<sup>2</sup> of the vector.
    */
    __forceinline float length2() const;
    
    //! Length.
    /*!
        Returns the geometric length of the vector.
    */
    __forceinline float length() const {
		ALIGN_SSE const float length_2 = length2();
		ALIGN_SSE float length_recip, l;
#ifndef NO_SSE
		fastrsqrtss(setSSE(length_2), length_recip);
		recipss(setSSE(length_recip), l);
#else
		l = sqrtf(length_2);
#endif
		return l;}
    
    //! Normalizes the vector and return its length.
    /*!
        Scales each component of the vector in order to get unit
        length without changing direction.
    
        \return The length of the vector prior to normalization.
    */
    __forceinline float unitize()
    {
        ALIGN_SSE float l;
#ifndef NO_SSE
		fastrsqrtss(setSSE(length2()), l);
#else
		l = 1.0f / sqrtf(length2());
#endif
        *this *= l;
        return l;
    }

	__forceinline float lengthRecip()
	{
		ALIGN_SSE float l;
#ifndef NO_SSE
		fastrsqrtss(setSSE(length2()), l);
#else
		l = 1.0f / sqrtf(length2());
#endif
		return l;
	}
    
    //! Normalize a vector and return a reference to it.
    /*!
        Scales each component of the vector in order to get unit
        length without changing direction.
    
        \return A reference to the vector.
    */
    __forceinline const Vector3 & normalize()
    {
		ALIGN_SSE float l;
#ifndef NO_SSE
		fastrsqrtss(setSSE(length2()), l);
#else
		l = 1.0f / sqrtf(length2());
#endif
        return (*this *= l);
    }
    
    //! Return a normalized copy of the vector.
    __forceinline Vector3 normalized() const
    {
		ALIGN_SSE float l;
#ifndef NO_SSE
		fastrsqrtss(setSSE(length2()), l);
#else
		l = 1.0f / sqrtf(length2());
#endif
        return( *this * l);
    }
    
    //! Return a rotated copy of the vector
    __forceinline Vector3 rotated(float theta, const Vector3 & w) const;
    
    //! Rotate this vector about another vector, w, by theta radians.
    const Vector3 & rotate(float theta, const Vector3 & w)
    {
	return *this = rotated(theta, w);
    }

	__forceinline const float average()
	{
		return (x+y+z) * 0.333333f;
	}

	__forceinline const float maxComp()
	{
		return std::max(x, std::max(y, z));
	}
};


//! Multiply a scalar by a Vec3.
__forceinline Vector3 operator*(float s, const Vector3& a)
{
	return Vector3(a.x * s, a.y * s, a.z * s);
}


//! The dot product of two Vec3s.
__forceinline float dot(const Vector3 & a, const Vector3 & b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


//! The cross product of two Vec3s.
__forceinline Vector3 cross(const Vector3 & a, const Vector3 & b)
{
    return Vector3(a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x);
}


__forceinline float Vector3::length2() const
{
    return dot(*this, *this);
}


//! Return a rotated copy of the vector
__forceinline Vector3 Vector3::rotated(float theta, const Vector3 & w) const
{
    float c = cosf(theta);
    float s = sinf(theta);

    Vector3 v0 = dot(*this, w) * w;
    Vector3 v1 = *this - v0;
    Vector3 v2 = cross(w, v1);

    return v0 + c * v1 + s * v2;
}


__forceinline std::ostream &
operator<<(std::ostream& out, const Vector3& a)
{
	return out << a.x << " " << a.y << " " << a.z ;
}

#endif // CSE168_VECTOR3_H_INCLUDED
