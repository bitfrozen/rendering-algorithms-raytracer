#ifndef CSE168_MATRIX4X4_H_INCLUDED
#define CSE168_MATRIX4X4_H_INCLUDED

#include <math.h>
#include <float.h>
#include <iostream>
#include "SSE.h"

#ifdef WIN32
#pragma warning(disable:4305) // disable useless warnings
#pragma warning(disable:4244)
#endif

#include "Vector3.h"
#include "Vector4.h"

ALIGN_SSE class Matrix4x4
{

public:
	union {float m1[4]; __m128 _m1; struct { float m11, m12, m13, m14; }; };
	union {float m2[4]; __m128 _m2; struct { float m21, m22, m23, m24; }; };
	union {float m3[4]; __m128 _m3; struct { float m31, m32, m33, m34; }; };
	union {float m4[4]; __m128 _m4; struct { float m41, m42, m43, m44; }; };

    // Implements a 4x4 matrix: m_ij - row-i and column-j entry

public:

    Matrix4x4();
	Matrix4x4(const Matrix4x4& M);
    Matrix4x4(const Vector4&, const Vector4&,
              const Vector4&, const Vector4&); // sets by columns!
	Matrix4x4(const __m128&, const __m128&,
			  const __m128&, const __m128&); // sets by rows!
    Matrix4x4(float, float, float, float,
              float, float, float, float,
              float, float, float, float,
              float, float, float, float); // sets by columns

    inline void setIdentity();     // set to the identity map
    inline void set(const Matrix4x4&); // set to the matrix.
    inline void set(const Vector4&, const Vector4&,
                    const Vector4&, const Vector4&);
    inline void set(float, float, float, float,
                    float, float, float, float,
                    float, float, float, float,
                    float, float, float, float);
    inline void setColumn1(float, float, float, float);
    inline void setColumn2(float, float, float, float);
    inline void setColumn3(float, float, float, float);
    inline void setColumn4(float, float, float, float);
    inline void setColumn1(const Vector4&);
    inline void setColumn2(const Vector4&);
    inline void setColumn3(const Vector4&);
    inline void setColumn4(const Vector4&);
    inline Vector4 column1() const;
    inline Vector4 column2() const;
    inline Vector4 column3() const;
    inline Vector4 column4() const;
    inline Matrix4x4& invert();

    inline void transpose();                    // Transposes it.
    inline Matrix4x4& operator+=(const Matrix4x4&);
    inline Matrix4x4& operator-=(const Matrix4x4&);
    inline Matrix4x4& operator*=(float);
    inline Matrix4x4& operator/=(float);
    inline Matrix4x4& operator*=(const Matrix4x4&);    // Matrix product

	const inline Vector4 multiplyAndDivideByW(const Vector4&) const;
	const inline Vector3 multiplyAndDivideByW(Vector3&) const;
#ifndef NO_SSE
	const inline Vector3 multiplyAndDivideByW(const __m128&) const;
#endif

	inline void translate(float x, float y, float z);
	inline void scale(float x, float y, float z);
	inline void rotate(float angle, float x, float y, float z);
};

inline std::ostream &
operator<<(std::ostream& out, const Matrix4x4& m)
{
    return out << m.column1() << std::endl
               << m.column2() << std::endl
               << m.column3() << std::endl
               << m.column4();
}

// binary operators
inline Matrix4x4 operator+(const Matrix4x4&, const Matrix4x4&);
inline Matrix4x4 operator-(const Matrix4x4&);
inline Matrix4x4 operator-(const Matrix4x4&, const Matrix4x4&);
inline Matrix4x4 operator*(const Matrix4x4&, float);
inline Matrix4x4 operator*(float, const Matrix4x4&);
inline Matrix4x4 operator/(const Matrix4x4&, float);
inline Matrix4x4 operator*(const Matrix4x4&, const Matrix4x4&);

inline Vector4 operator*(const Matrix4x4&, const Vector4&);
inline Vector3 operator*(const Matrix4x4&, const Vector3&);



inline Matrix4x4::Matrix4x4() {setIdentity();}

inline Matrix4x4::Matrix4x4(const __m128& _M1, const __m128& _M2, const __m128& _M3, const __m128& _M4)
{
	_m1 = _M1;
	_m2 = _M2;
	_m3 = _M3;
	_m4 = _M4;
}

inline Matrix4x4::Matrix4x4(const Matrix4x4& M) // set to the matrix.
{
	m11 = M.m11;
	m12 = M.m12;
	m13 = M.m13;
	m14 = M.m14;
	m21 = M.m21;
	m22 = M.m22;
	m23 = M.m23;
	m24 = M.m24;
	m31 = M.m31;
	m32 = M.m32;
	m33 = M.m33;
	m34 = M.m34;
	m41 = M.m41;
	m42 = M.m42;
	m43 = M.m43;
	m44 = M.m44;
}

inline Matrix4x4::Matrix4x4(const Vector4& u, const Vector4& v,
                            const Vector4& s, const Vector4& t)
{
    m11 = u.x;      // Column 1
    m21 = u.y;
    m31 = u.z;
    m41 = u.w;
    m12 = v.x;      // Column 2
    m22 = v.y;
    m32 = v.z;
    m42 = v.w;
    m13 = s.x;      // Column 3
    m23 = s.y;
    m33 = s.z;
    m43 = s.w;
    m14 = t.x;      // Column 4
    m24 = t.y;
    m34 = t.z;
    m44 = t.w;
}

inline Matrix4x4::Matrix4x4(float a11, float a12, float a13, float a14,
                            float a21, float a22, float a23, float a24,
                            float a31, float a32, float a33, float a34,
                            float a41, float a42, float a43, float a44)
// Values specified in row order!!!
{
    m11 = a11;      // Row 1
    m12 = a12;
    m13 = a13;
    m14 = a14;
    m21 = a21;      // Row 2
    m22 = a22;
    m23 = a23;
    m24 = a24;
    m31 = a31;      // Row 3
    m32 = a32;
    m33 = a33;
    m34 = a34;
    m41 = a41;      // Row 4
    m42 = a42;
    m43 = a43;
    m44 = a44;
}

inline void
Matrix4x4::setIdentity()
{
    m11 = m22 = m33 = m44 = 1.0;
    m12 = m13 = m14 = m21 = m23 = m24 = m31 = m32 = m34 = m41= m42 = m43 = 0.0;
}

inline void
Matrix4x4::set(const Vector4& u, const Vector4& v,
               const Vector4& s, const Vector4& t)
{
    m11 = u.x;      // Column 1
    m21 = u.y;
    m31 = u.z;
    m41 = u.w;
    m12 = v.x;      // Column 2
    m22 = v.y;
    m32 = v.z;
    m42 = v.w;
    m13 = s.x;      // Column 3
    m23 = s.y;
    m33 = s.z;
    m43 = s.w;
    m14 = t.x;      // Column 4
    m24 = t.y;
    m34 = t.z;
    m44 = t.w;
}

inline void
Matrix4x4::set(float a11, float a12, float a13, float a14,
               float a21, float a22, float a23, float a24,
               float a31, float a32, float a33, float a34,
               float a41, float a42, float a43, float a44)
// Values specified in row order!!!
{
    m11 = a11;      // Row 1
    m12 = a12;
    m13 = a13;
    m14 = a14;
    m21 = a21;      // Row 2
    m22 = a22;
    m23 = a23;
    m24 = a24;
    m31 = a31;      // Row 3
    m32 = a32;
    m33 = a33;
    m34 = a34;
    m41 = a41;      // Row 4
    m42 = a42;
    m43 = a43;
    m44 = a44;
}

inline void
Matrix4x4::set(const Matrix4x4& M) // set to the matrix.
{
    m11 = M.m11;
    m12 = M.m12;
    m13 = M.m13;
    m14 = M.m14;
    m21 = M.m21;
    m22 = M.m22;
    m23 = M.m23;
    m24 = M.m24;
    m31 = M.m31;
    m32 = M.m32;
    m33 = M.m33;
    m34 = M.m34;
    m41 = M.m41;
    m42 = M.m42;
    m43 = M.m43;
    m44 = M.m44;
}

inline void
Matrix4x4::setColumn1(float x, float y, float z, float w)
{
    m11 = x; m21 = y; m31= z; m41 = w;
}

inline void
Matrix4x4::setColumn2(float x, float y, float z, float w)
{
    m12 = x; m22 = y; m32= z; m42 = w;
}

inline void
Matrix4x4::setColumn3(float x, float y, float z, float w)
{
    m13 = x; m23 = y; m33= z; m43 = w;
}

inline void
Matrix4x4::setColumn4(float x, float y, float z, float w)
{
    m14 = x; m24 = y; m34= z; m44 = w;
}

inline void
Matrix4x4::setColumn1(const Vector4& u)
{
    m11 = u.x; m21 = u.y; m31 = u.z; m41 = u.w;
}

inline void
Matrix4x4::setColumn2(const Vector4& u)
{
    m12 = u.x; m22 = u.y; m32 = u.z; m42 = u.w;
}

inline void
Matrix4x4::setColumn3(const Vector4& u)
{
    m13 = u.x; m23 = u.y; m33 = u.z; m43 = u.w;
}

inline void
Matrix4x4::setColumn4(const Vector4& u)
{
    m14 = u.x; m24 = u.y; m34 = u.z; m44 = u.w;
}

Vector4
Matrix4x4::column1() const
{
    return Vector4(m11, m21, m31, m41);
}

Vector4
Matrix4x4::column2() const
{
    return Vector4(m12, m22, m32, m42);
}

Vector4
Matrix4x4::column3() const
{
    return Vector4(m13, m23, m33, m43);
}

Vector4
Matrix4x4::column4() const
{
    return Vector4(m14, m24, m34, m44);
}

inline void
Matrix4x4::transpose()  // Transposes it.
{
    register float temp;
    temp = m12;
    m12 = m21;
    m21=temp;
    temp = m13;
    m13 = m31;
    m31 = temp;
    temp = m14;
    m14 = m41;
    m41 = temp;
    temp = m23;
    m23 = m32;
    m32 = temp;
    temp = m24;
    m24 = m42;
    m42 = temp;
    temp = m34;
    m34 = m43;
    m43 = temp;
}

Matrix4x4&
Matrix4x4::invert()          // Converts into inverse.
{
    float Tbt34C12 = m31*m42-m32*m41;       // 2x2 subdeterminants
    float Tbt34C13 = m31*m43-m33*m41;
    float Tbt34C14 = m31*m44-m34*m41;
    float Tbt34C23 = m32*m43-m33*m42;
    float Tbt34C24 = m32*m44-m34*m42;
    float Tbt34C34 = m33*m44-m34*m43;
    float Tbt24C12 = m21*m42-m22*m41;       // 2x2 subdeterminants
    float Tbt24C13 = m21*m43-m23*m41;
    float Tbt24C14 = m21*m44-m24*m41;
    float Tbt24C23 = m22*m43-m23*m42;
    float Tbt24C24 = m22*m44-m24*m42;
    float Tbt24C34 = m23*m44-m24*m43;
    float Tbt23C12 = m21*m32-m22*m31;       // 2x2 subdeterminants
    float Tbt23C13 = m21*m33-m23*m31;
    float Tbt23C14 = m21*m34-m24*m31;
    float Tbt23C23 = m22*m33-m23*m32;
    float Tbt23C24 = m22*m34-m24*m32;
    float Tbt23C34 = m23*m34-m24*m33;

    float sd11 = m22*Tbt34C34 - m23*Tbt34C24 + m24*Tbt34C23;    // 3x3 subdeterminants
    float sd12 = m21*Tbt34C34 - m23*Tbt34C14 + m24*Tbt34C13;
    float sd13 = m21*Tbt34C24 - m22*Tbt34C14 + m24*Tbt34C12;
    float sd14 = m21*Tbt34C23 - m22*Tbt34C13 + m23*Tbt34C12;
    float sd21 = m12*Tbt34C34 - m13*Tbt34C24 + m14*Tbt34C23;    // 3x3 subdeterminants
    float sd22 = m11*Tbt34C34 - m13*Tbt34C14 + m14*Tbt34C13;
    float sd23 = m11*Tbt34C24 - m12*Tbt34C14 + m14*Tbt34C12;
    float sd24 = m11*Tbt34C23 - m12*Tbt34C13 + m13*Tbt34C12;
    float sd31 = m12*Tbt24C34 - m13*Tbt24C24 + m14*Tbt24C23;    // 3x3 subdeterminants
    float sd32 = m11*Tbt24C34 - m13*Tbt24C14 + m14*Tbt24C13;
    float sd33 = m11*Tbt24C24 - m12*Tbt24C14 + m14*Tbt24C12;
    float sd34 = m11*Tbt24C23 - m12*Tbt24C13 + m13*Tbt24C12;
    float sd41 = m12*Tbt23C34 - m13*Tbt23C24 + m14*Tbt23C23;    // 3x3 subdeterminants
    float sd42 = m11*Tbt23C34 - m13*Tbt23C14 + m14*Tbt23C13;
    float sd43 = m11*Tbt23C24 - m12*Tbt23C14 + m14*Tbt23C12;
    float sd44 = m11*Tbt23C23 - m12*Tbt23C13 + m13*Tbt23C12;

    register float detInv = 1.0/(m11*sd11 - m12*sd12 + m13*sd13 - m14*sd14);

    m11 = sd11*detInv;
    m12 = -sd21*detInv;
    m13 = sd31*detInv;
    m14 = -sd41*detInv;
    m21 = -sd12*detInv;
    m22 = sd22*detInv;
    m23 = -sd32*detInv;
    m24 = sd42*detInv;
    m31 = sd13*detInv;
    m32 = -sd23*detInv;
    m33 = sd33*detInv;
    m34 = -sd43*detInv;
    m41 = -sd14*detInv;
    m42 = sd24*detInv;
    m43 = -sd34*detInv;
    m44 = sd44*detInv;

    return *this;
}

inline Matrix4x4&
Matrix4x4::operator+=(const Matrix4x4& B)
{
#ifndef NO_SSE
    _m1 = addps(_m1, B._m1);
	_m2 = addps(_m2, B._m2);
	_m3 = addps(_m3, B._m3);
	_m4 = addps(_m4, B._m4);
#else
	m11 += B.m11;
    m12 += B.m12;
    m13 += B.m13;
    m14 += B.m14;
    m21 += B.m21;
    m22 += B.m22;
    m23 += B.m23;
    m24 += B.m24;
    m31 += B.m31;
    m32 += B.m32;
    m33 += B.m33;
    m34 += B.m34;
    m41 += B.m41;
    m42 += B.m42;
    m43 += B.m43;
    m44 += B.m44;
#endif
    return *this;
}

inline Matrix4x4&
Matrix4x4::operator-=(const Matrix4x4& B)
{
#ifndef NO_SSE
	_m1 = subps(_m1, B._m1);
	_m2 = subps(_m2, B._m2);
	_m3 = subps(_m3, B._m3);
	_m4 = subps(_m4, B._m4);
#else
    m11 -= B.m11;
    m12 -= B.m12;
    m13 -= B.m13;
    m14 -= B.m14;
    m21 -= B.m21;
    m22 -= B.m22;
    m23 -= B.m23;
    m24 -= B.m24;
    m31 -= B.m31;
    m32 -= B.m32;
    m33 -= B.m33;
    m34 -= B.m34;
    m41 -= B.m41;
    m42 -= B.m42;
    m43 -= B.m43;
    m44 -= B.m44;
#endif
    return(*this);
}

inline Matrix4x4
operator+(const Matrix4x4& A, const Matrix4x4& B)
{
#ifndef NO_SSE
	return Matrix4x4(addps(A._m1, B._m1), addps(A._m2, B._m2), addps(A._m3, B._m3), addps(A._m4, B._m4));
#else
    return Matrix4x4(A.m11+B.m11, A.m21+B.m21, A.m31+B.m31, A.m41+B.m41,
                     A.m12+B.m12, A.m22+B.m22, A.m32+B.m32, A.m42+B.m42,
                     A.m13+B.m13, A.m23+B.m23, A.m33+B.m33, A.m43+B.m43,
                     A.m14+B.m14, A.m24+B.m24, A.m34+B.m34, A.m44+B.m44);
#endif
}

inline Matrix4x4
operator-(const Matrix4x4& A)
{
    return Matrix4x4(-A.m11, -A.m21, -A.m31, -A.m41,
                     -A.m12, -A.m22, -A.m32, -A.m42,
                     -A.m13, -A.m23, -A.m33, -A.m43,
                     -A.m14, -A.m24, -A.m34, -A.m44);
}

inline Matrix4x4
operator-(const Matrix4x4& A, const Matrix4x4& B)
{
#ifndef NO_SSE
	return Matrix4x4(subps(A._m1, B._m1), subps(A._m2, B._m2), subps(A._m3, B._m3), subps(A._m4, B._m4));
#else
    return Matrix4x4(A.m11-B.m11, A.m21-B.m21, A.m31-B.m31, A.m41-B.m41,
                     A.m12-B.m12, A.m22-B.m22, A.m32-B.m32, A.m42-B.m42,
                     A.m13-B.m13, A.m23-B.m23, A.m33-B.m33, A.m43-B.m43,
                     A.m14-B.m14, A.m24-B.m24, A.m34-B.m34, A.m44-B.m44);
#endif
}

inline Matrix4x4&
Matrix4x4::operator*=(float b)
{
#ifndef NO_SSE
	const __m128 f = setSSE(b);
	_m1 = mulps(_m1, f);
	_m2 = mulps(_m2, f);
	_m3 = mulps(_m3, f);
	_m4 = mulps(_m4, f);
#else
    m11 *= b;
    m12 *= b;
    m13 *= b;
    m14 *= b;
    m21 *= b;
    m22 *= b;
    m23 *= b;
    m24 *= b;
    m31 *= b;
    m32 *= b;
    m33 *= b;
    m34 *= b;
    m41 *= b;
    m42 *= b;
    m43 *= b;
    m44 *= b;
#endif
    return *this;
}

inline Matrix4x4&
Matrix4x4::operator*=(const Matrix4x4& B)    // Matrix product
{
#ifndef NO_SSE
	const __m128 _b1 = _mm_set_ps(B.m41, B.m31, B.m21, B.m11);
	const __m128 _b2 = _mm_set_ps(B.m42, B.m32, B.m22, B.m12);
	const __m128 _b3 = _mm_set_ps(B.m43, B.m33, B.m23, B.m13);
	const __m128 _b4 = _mm_set_ps(B.m44, B.m34, B.m24, B.m14);

	_m1 = shuffleps(_mm_movelh_ps(dotps(_m1, _b1, 0xFF), dotps(_m1, _b2, 0xFF)), _mm_movelh_ps(dotps(_m1, _b3, 0xFF), dotps(_m1, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	_m2 = shuffleps(_mm_movelh_ps(dotps(_m2, _b1, 0xFF), dotps(_m2, _b2, 0xFF)), _mm_movelh_ps(dotps(_m2, _b3, 0xFF), dotps(_m2, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	_m3 = shuffleps(_mm_movelh_ps(dotps(_m3, _b1, 0xFF), dotps(_m3, _b2, 0xFF)), _mm_movelh_ps(dotps(_m3, _b3, 0xFF), dotps(_m3, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	_m4 = shuffleps(_mm_movelh_ps(dotps(_m4, _b1, 0xFF), dotps(_m4, _b2, 0xFF)), _mm_movelh_ps(dotps(_m4, _b3, 0xFF), dotps(_m4, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
#else
    float t1, t2, t3;       // temporary values
    t1 =  m11*B.m11 + m12*B.m21 + m13*B.m31 + m14*B.m41;
    t2 =  m11*B.m12 + m12*B.m22 + m13*B.m32 + m14*B.m42;
    t3 =  m11*B.m13 + m12*B.m23 + m13*B.m33 + m14*B.m43;
    m14 = m11*B.m14 + m12*B.m24 + m13*B.m34 + m14*B.m44;
    m11 = t1;
    m12 = t2;
    m13 = t3;

    t1 =  m21*B.m11 + m22*B.m21 + m23*B.m31 + m24*B.m41;
    t2 =  m21*B.m12 + m22*B.m22 + m23*B.m32 + m24*B.m42;
    t3 =  m21*B.m13 + m22*B.m23 + m23*B.m33 + m24*B.m43;
    m24 = m21*B.m14 + m22*B.m24 + m23*B.m34 + m24*B.m44;
    m21 = t1;
    m22 = t2;
    m23 = t3;

    t1 =  m31*B.m11 + m32*B.m21 + m33*B.m31 + m34*B.m41;
    t2 =  m31*B.m12 + m32*B.m22 + m33*B.m32 + m34*B.m42;
    t3 =  m31*B.m13 + m32*B.m23 + m33*B.m33 + m34*B.m43;
    m34 = m31*B.m14 + m32*B.m24 + m33*B.m34 + m34*B.m44;
    m31 = t1;
    m32 = t2;
    m33 = t3;

    t1 =  m41*B.m11 + m42*B.m21 + m43*B.m31 + m44*B.m41;
    t2 =  m41*B.m12 + m42*B.m22 + m43*B.m32 + m44*B.m42;
    t3 =  m41*B.m13 + m42*B.m23 + m43*B.m33 + m44*B.m43;
    m44 = m41*B.m14 + m42*B.m24 + m43*B.m34 + m44*B.m44;
    m41 = t1;
    m42 = t2;
    m43 = t3;
#endif
    return *this;
}

inline Matrix4x4
operator*(const Matrix4x4& A, const Matrix4x4& B) // Matrix product
{
	Matrix4x4 R;
#ifndef NO_SSE
	const __m128 _b1 = _mm_set_ps(B.m41, B.m31, B.m21, B.m11);
	const __m128 _b2 = _mm_set_ps(B.m42, B.m32, B.m22, B.m12);
	const __m128 _b3 = _mm_set_ps(B.m43, B.m33, B.m23, B.m13);
	const __m128 _b4 = _mm_set_ps(B.m44, B.m34, B.m24, B.m14);

	R._m1 = shuffleps(_mm_movelh_ps(dotps(A._m1, _b1, 0xFF), dotps(A._m1, _b2, 0xFF)), _mm_movelh_ps(dotps(A._m1, _b3, 0xFF), dotps(A._m1, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	R._m2 = shuffleps(_mm_movelh_ps(dotps(A._m2, _b1, 0xFF), dotps(A._m2, _b2, 0xFF)), _mm_movelh_ps(dotps(A._m2, _b3, 0xFF), dotps(A._m2, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	R._m3 = shuffleps(_mm_movelh_ps(dotps(A._m3, _b1, 0xFF), dotps(A._m3, _b2, 0xFF)), _mm_movelh_ps(dotps(A._m3, _b3, 0xFF), dotps(A._m3, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
	R._m4 = shuffleps(_mm_movelh_ps(dotps(A._m4, _b1, 0xFF), dotps(A._m4, _b2, 0xFF)), _mm_movelh_ps(dotps(A._m4, _b3, 0xFF), dotps(A._m4, _b4, 0xFF)), _MM_SHUFFLE(2,0,2,0));
#else
    float t1, t2, t3;       // temporary values
    t1 =  A.m11*B.m11 + A.m12*B.m21 + A.m13*B.m31 + A.m14*B.m41;
    t2 =  A.m11*B.m12 + A.m12*B.m22 + A.m13*B.m32 + A.m14*B.m42;
    t3 =  A.m11*B.m13 + A.m12*B.m23 + A.m13*B.m33 + A.m14*B.m43;
    R.m14 = A.m11*B.m14 + A.m12*B.m24 + A.m13*B.m34 + A.m14*B.m44;
    R.m11 = t1;
    R.m12 = t2;
    R.m13 = t3;

    t1 =  A.m21*B.m11 + A.m22*B.m21 + A.m23*B.m31 + A.m24*B.m41;
    t2 =  A.m21*B.m12 + A.m22*B.m22 + A.m23*B.m32 + A.m24*B.m42;
    t3 =  A.m21*B.m13 + A.m22*B.m23 + A.m23*B.m33 + A.m24*B.m43;
    R.m24 = A.m21*B.m14 + A.m22*B.m24 + A.m23*B.m34 + A.m24*B.m44;
    R.m21 = t1;
    R.m22 = t2;
    R.m23 = t3;

    t1 =  A.m31*B.m11 + A.m32*B.m21 + A.m33*B.m31 + A.m34*B.m41;
    t2 =  A.m31*B.m12 + A.m32*B.m22 + A.m33*B.m32 + A.m34*B.m42;
    t3 =  A.m31*B.m13 + A.m32*B.m23 + A.m33*B.m33 + A.m34*B.m43;
    R.m34 = A.m31*B.m14 + A.m32*B.m24 + A.m33*B.m34 + A.m34*B.m44;
    R.m31 = t1;
    R.m32 = t2;
    R.m33 = t3;

    t1 =  A.m41*B.m11 + A.m42*B.m21 + A.m43*B.m31 + A.m44*B.m41;
    t2 =  A.m41*B.m12 + A.m42*B.m22 + A.m43*B.m32 + A.m44*B.m42;
    t3 =  A.m41*B.m13 + A.m42*B.m23 + A.m43*B.m33 + A.m44*B.m43;
    R.m44 = A.m41*B.m14 + A.m42*B.m24 + A.m43*B.m34 + A.m44*B.m44;
    R.m41 = t1;
    R.m42 = t2;
    R.m43 = t3;
#endif
    return R;
}

inline Matrix4x4
operator*(const Matrix4x4& A, float b)
{
#ifndef NO_SSE
	const __m128 f = setSSE(b);
	return(Matrix4x4(mulps(A._m1, f), mulps(A._m2, f), mulps(A._m3, f), mulps(A._m4, f)));
#else
    return(Matrix4x4(A.m11*b, A.m21*b, A.m31*b, A.m41*b,
                     A.m12*b, A.m22*b, A.m32*b, A.m42*b,
                     A.m13*b, A.m23*b, A.m33*b, A.m43*b,
                     A.m14*b, A.m24*b, A.m34*b, A.m44*b));
#endif
}

inline Matrix4x4
operator*(float b, const Matrix4x4& A)
{
#ifndef NO_SSE
	const __m128 f = setSSE(b);
	return(Matrix4x4(mulps(A._m1, f), mulps(A._m2, f), mulps(A._m3, f), mulps(A._m4, f)));
#else
    return(Matrix4x4(A.m11*b, A.m21*b, A.m31*b, A.m41*b,
                     A.m12*b, A.m22*b, A.m32*b, A.m42*b,
                     A.m13*b, A.m23*b, A.m33*b, A.m43*b,
                     A.m14*b, A.m24*b, A.m34*b, A.m44*b));
#endif
}

inline Matrix4x4
operator/(const Matrix4x4& A, float b)
{
    register float bInv = 1.0/b;
    return (A*bInv);
}

inline Matrix4x4&
Matrix4x4::operator/=(float b)
{
    register float bInv = 1.0/b;
    return (*this *= bInv);
}

inline Vector4
operator*(const Matrix4x4& A, const Vector4& u)
{
#ifndef NO_SSE
	return Vector4(shuffleps(_mm_movelh_ps(dotps(A._m1, u._v, 0xFF), dotps(A._m2, u._v, 0xFF)), _mm_movelh_ps(dotps(A._m3, u._v, 0xFF), dotps(A._m4, u._v, 0xFF)), _MM_SHUFFLE(2,0,2,0)));
#else
    return Vector4(A.m11*u.x + A.m12*u.y + A.m13*u.z + A.m14*u.w,
                   A.m21*u.x + A.m22*u.y + A.m23*u.z + A.m24*u.w,
                   A.m31*u.x + A.m32*u.y + A.m33*u.z + A.m34*u.w,
                   A.m41*u.x + A.m42*u.y + A.m43*u.z + A.m44*u.w);
#endif
}

inline Vector3
operator*(const Matrix4x4& A, const Vector3& u)
{
/*#ifndef NO_SSE
	return Vector3(shuffleps(_mm_movelh_ps(dotps(A._m1, u._v, 0xFF), dotps(A._m2, u._v, 0xFF)), _mm_movelh_ps(dotps(A._m3, u._v, 0xFF), setZero), _MM_SHUFFLE(2,0,2,0)));
#else*/
    return Vector3(A.m11*u.x + A.m12*u.y + A.m13*u.z,
                   A.m21*u.x + A.m22*u.y + A.m23*u.z,
                   A.m31*u.x + A.m32*u.y + A.m33*u.z);
    // note that this ignores the fourth row in the matrix!
//#endif
}

#ifndef NO_SSE
inline Vector3
operator*(const Matrix4x4& A, const __m128& u)
{
	return Vector3(shuffleps(_mm_movelh_ps(dotps(A._m1, u, 0xFF), dotps(A._m2, u, 0xFF)), _mm_movelh_ps(dotps(A._m3, u, 0xFF), setZero), _MM_SHUFFLE(2,0,2,0)));
}
#endif

const inline Vector4 Matrix4x4::multiplyAndDivideByW(const Vector4& u) const
{
#ifndef NO_SSE
	const __m128 wRecip = recipps(dotps(_m4, u._v, 0xFF));
	return Vector4(mulps(wRecip, shuffleps(_mm_movelh_ps(dotps(_m1, u._v, 0xFF), dotps(_m2, u._v, 0xFF)), _mm_movelh_ps(dotps(_m3, u._v, 0xFF), _onesps), _MM_SHUFFLE(2,0,2,0))));
#else
	const float wRecip = 1.0f / (A.m41*u.x + A.m42*u.y + A.m43*u.z + A.m44*u.w);
	return Vector4((A.m11*u.x + A.m12*u.y + A.m13*u.z + A.m14*u.w) * wRecip,
				   (A.m21*u.x + A.m22*u.y + A.m23*u.z + A.m24*u.w) * wRecip,
				   (A.m31*u.x + A.m32*u.y + A.m33*u.z + A.m34*u.w) * wRecip,
				   1.0f);
#endif
}

const inline Vector3 Matrix4x4::multiplyAndDivideByW(Vector3& u) const
{
#ifndef NO_SSE
	u.__dummy = 1.0f;
	const __m128 wRecip = recipps(dotps(_m4, u._v, 0xFF));
	return Vector3(mulps(wRecip, shuffleps(_mm_movelh_ps(dotps(_m1, u._v, 0xFF), dotps(_m2, u._v, 0xFF)), _mm_movelh_ps(dotps(_m3, u._v, 0xFF), setZero), _MM_SHUFFLE(2,0,2,0))));
#else
	u.__dummy = 1.0f;
	const float wRecip = 1.0f / (A.m41*u.x + A.m42*u.y + A.m43*u.z + A.m44);
	return Vector3((A.m11*u.x + A.m12*u.y + A.m13*u.z + A.m14) * wRecip,
				   (A.m21*u.x + A.m22*u.y + A.m23*u.z + A.m24) * wRecip,
				   (A.m31*u.x + A.m32*u.y + A.m33*u.z + A.m34) * wRecip);
#endif
}

#ifndef NO_SSE
const inline Vector3 Matrix4x4::multiplyAndDivideByW(const __m128& u) const
{
	const __m128 wRecip = recipps(dotps(_m4, u, 0xFF));
	return Vector3(mulps(wRecip, shuffleps(_mm_movelh_ps(dotps(_m1, u, 0xFF), dotps(_m2, u, 0xFF)), _mm_movelh_ps(dotps(_m3, u, 0xFF), setZero), _MM_SHUFFLE(2,0,2,0))));
}
#endif

inline void Matrix4x4::translate(float x, float y, float z)
{
	setColumn4(Vector4(x, y, z, 0) + Vector4(m14, m24, m34, m44));
}


inline void Matrix4x4::scale(float x, float y, float z)
{
	m11 *= x;
	m22 *= y;
	m33 *= z;
}

// angle is in degrees
inline void Matrix4x4::rotate(float angle, float x, float y, float z)
{
	float rad = angle*(PI/180.);

	float x2 = x*x;
	float y2 = y*y;
	float z2 = z*z;
	float c = cos(rad);
	float cinv = 1-c;
	float s = sin(rad);
	float xy = x*y;
	float xz = x*z;
	float yz = y*z;
	float xs = x*s;
	float ys = y*s;
	float zs = z*s;
	float xzcinv = xz*cinv;
	float xycinv = xy*cinv;
	float yzcinv = yz*cinv;

	set(x2 + c*(1-x2), xy*cinv+zs, xzcinv - ys, 0,
		xycinv - zs, y2 + c*(1-y2), yzcinv + xs, 0,
		xzcinv + ys, yzcinv - xs, z2 + c*(1-z2), 0,
		0, 0, 0, 1);
}

#endif // CSE168_MATRIX4X4_H_INCLUDED