
#include "RawImage.h"

RawImage::RawImage(int w, int h, float* data, ImageType t) {
	m_width = w;
	m_height = h;
	m_rawData = data;
	m_imageType = t;
}