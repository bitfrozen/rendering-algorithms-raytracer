#include "Texture.h"
#include "Miro.h"

Texture::Texture() 
{
}

Texture::~Texture() 
{
}

const float Texture::getLookupAlpha(float u, float v) const
{
	//force to be between 0 and 1
	u = u - float(int(u));
	v = v - float(int(v));
	if (u < 0.0f)
		u = u + 1.0f;
	if (v < 0.0f)
		v = v + 1.0f;

	//textures start with v = 0 at the top, so reverse our v:
	v = 1.0f - v;

	//use bilinear lookup
	float px = u*(m_image->m_width);
	float py = v*(m_image->m_height);
	float x1 = floor(px);
	float x2 = x1 + 1.0f;
	float dx = px - x1;
	float y1 = floor(py);
	float y2 = y1 + 1.0f;
	float dy = py - y1;

	//interpolate in the x direction:
	Vector4 q1 = getPixel(x1, y1)*(1.0f - dx) + getPixel(x2, y1)*dx;
	Vector4 q2 = getPixel(x1, y2)*(1.0f - dx) + getPixel(x2, y2)*dx;

	//interpolate those two points
	return q1.w*(1.0f - dy) + q2.w*dy;
}

const Vector4 Texture::getLookup(float u, float v) const
{
	//force to be between 0 and 1
	u = u - float(int(u));
	v = v - float(int(v));
	if (u < 0.0f)
		u = u + 1.0f;
	if (v < 0.0f)
		v = v + 1.0f;

	//textures start with v = 0 at the top, so reverse our v:
	v = 1.0f - v;

    //use bilinear lookup
	float px = u*(m_image->m_width);
	float py = v*(m_image->m_height);
	float x1 = floor(px);
    float x2 = x1 + 1.0f;
	float dx = px - x1;
	float y1 = floor(py);
	float y2 = y1 + 1.0f;
	float dy = py - y1;

	//interpolate in the x direction:
	Vector4 q1 = getPixel(x1, y1)*(1.0f - dx) + getPixel(x2, y1)*dx;
	Vector4 q2 = getPixel(x1, y2)*(1.0f - dx) + getPixel(x2, y2)*dx;

	//interpolate those two points
	return q1*(1.0f - dy) + q2*dy;
}

const Vector3 Texture::getLookup3(float u, float v) const
{
	Vector4 temp = getLookup(u,v);
	return Vector3(temp[0], temp[1], temp[2]);
}

const Vector3 Texture::getLookupXYZ3(float x, float y, float z) const 
{
	float theta = atan2(z, x) + PI;
	float phi = acos(y);
	float u = theta * 0.5 * _1_PI;
	float v = 1.0 - (phi * _1_PI);

	return getLookup3(u,v);
}

const Vector3 Texture::getLookupXYZ3(const Vector3& direction) const 
{
	float theta = atan2(direction.z, direction.x) + PI;
	float phi = acos(direction.y);
	float u = theta * 0.5 * _1_PI;
	float v = 1.0 - (phi * _1_PI);

	return getLookup3(u,v);
}

const Vector4 Texture::getPixel(int x, int y) const
{
	//clamp
	//if (x > m_image->m_width) x = m_image->m_width;
	//if (y > m_image->m_height) y = m_image->m_height;

	//tile
	x = x % m_image->m_width;
	y = y % m_image->m_height;
	
	if (m_image->m_imageType == GRAYSCALE) {
		float g = m_image->m_rawData[y*m_image->m_width + x];
		return Vector4(g, g, g, 1.0f);
	} else if (m_image->m_imageType == RGB) {
		int base = y*m_image->m_width*3 + x*3;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], 1.0f);
	} else if (m_image->m_imageType == RGBA) {
		int base = y*m_image->m_width*4 + x*4;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], m_image->m_rawData[base + 3]);
	} else if (m_image->m_imageType == HDR) {
		int base = y*m_image->m_width*3 + x*3;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], m_image->m_rawData[base + 3]);
	} else {
		return Vector4(0.0f, 0.0f, 0.0f, 1.0f); 
	}
}
