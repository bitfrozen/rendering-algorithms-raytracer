
#include "Texture.h"


Texture::Texture() {

}

Texture::~Texture() {
}

Vector4 Texture::getLookup(float u, float v) {
    //use bilinear lookup
	float px = u*m_image->m_width;
	float py = v*m_image->m_height;
	float x1 = floor(px);
    float x2 = x1 + 1;
	float dx = px - x1;
	float y1 = floor(py);
	float y2 = y1 + 1;
	float dy = py - y1;

	//interpolate in the x direction:
	Vector4 q1 = getPixel(x1, y1)*(1 - dx) + getPixel(x2, y1)*dx;
	Vector4 q2 = getPixel(x1, y2)*(1 - dx) + getPixel(x2, y2)*dx;

	//interpolate those two points
	return q1*(1 - dy) + q2*dy;
}

Vector4 Texture::getPixel(int x, int y) {
	if (m_image->m_imageType == GRAYSCALE) {
		float g = m_image->m_rawData[x*m_image->m_height + y];
		return Vector4(g, g, g, 1.0f);
	} else if (m_image->m_imageType == RGB) {
		int base = x*m_image->m_height*3 + y*3;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], 1.0f);
	} else if (m_image->m_imageType == RGBA) {
		int base = x*m_image->m_height*4 + y*4;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], m_image->m_rawData[base + 3]);
	} else if (m_image->m_imageType == HDR) {
        //TODO: How do I access HDR????
		int base = x*m_image->m_height*4 + y*4;
		return Vector4(m_image->m_rawData[base], m_image->m_rawData[base + 1], m_image->m_rawData[base + 2], m_image->m_rawData[base + 3]);
	} else {
		return Vector4(0.0f, 0.0f, 0.0f, 1.0f); 
	}
}
