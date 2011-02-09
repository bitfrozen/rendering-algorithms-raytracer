#ifndef CSE168_TEXTURE_INCLUDED
#define CSE168_TEXTURE_INCLUDED

#include "Vector4.h"
#include "Vector3.h"
#include "RawImage.h"


class Texture {

public:
	Texture();
	Texture(RawImage* image) { m_image = image; }
	~Texture();
	virtual Vector4 getLookup(float u, float v);
	virtual Vector3 getLookup3(float u, float v);
	/* Lookup for environment map. assumes texture is latitue-longitude map. */
	virtual Vector3 getLookupXYZ3(float x, float y, float z);

private:
	Vector4 getPixel(int x, int y);
	//pixels stored in column major format
	RawImage* m_image;
};
#endif