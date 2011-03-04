#ifndef CSE168_PROXYMATRIX_H_INCLUDED
#define CSE168_PROXYMATRIX_H_INCLUDED

#include "Matrix4x4.h"

ALIGN_SSE class ProxyMatrix 
{
public:
	ALIGN_SSE Matrix4x4 m_transform;
	ALIGN_SSE Matrix4x4 m_inverse;
	ALIGN_SSE Matrix4x4 m_invTranspose;
	ProxyMatrix(const Matrix4x4& M);
	ProxyMatrix();
	void set(const Matrix4x4& M);
};

#endif