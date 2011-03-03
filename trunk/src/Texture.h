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
	const virtual float getWidth() {return m_image->m_width;}
	const virtual float getHeight() {return m_image->m_height;}
	const virtual Vector4 getLookup(float u, float v) const;
	const virtual float   getLookupAlpha(float u, float v) const;
	const virtual Vector3 getLookup3(float u, float v) const;
	/* Lookup for environment map. assumes texture is latitude-longitude map. */
	const virtual Vector3 getLookupXYZ3(float x, float y, float z) const;
	const virtual Vector3 getLookupXYZ3(const Vector3& direction) const;

private:
	const Vector4 getPixel(int x, int y) const;
	//pixels stored in column major format
	RawImage* m_image;	
};
#endif