#include "ProxyMatrix.h"

ProxyMatrix::ProxyMatrix(const Matrix4x4& M)
{
	m_transform = M;
	m_inverse = M; m_inverse.invert();
	m_invTranspose = M; m_invTranspose.invert(); m_invTranspose.transpose();
}

ProxyMatrix::ProxyMatrix()
{
	m_transform = Matrix4x4();
	m_inverse = m_transform;
	m_invTranspose = m_transform;
}

void ProxyMatrix::set(const Matrix4x4& M)
{
	m_transform = M;
	m_inverse = M; m_inverse.invert();
	m_invTranspose = M; m_invTranspose.invert(); m_invTranspose.transpose();
}