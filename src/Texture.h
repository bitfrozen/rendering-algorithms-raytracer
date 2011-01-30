#ifndef CSE168_TEXTURE_INCLUDED
#define CSE168_TEXTURE_INCLUDED

#include "Vector4.h"
#include "RawImage.h"

class Texture {

public:
	Texture();
	Texture(RawImage* image) { m_image = image; }
	~Texture();
	Vector4 getLookup(float u, float v);

private:
	Vector4 getPixel(int x, int y);
	//pixels stored in column major format
	RawImage* m_image;
};
#endif