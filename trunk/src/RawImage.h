#ifndef CSE168_RAW_IMAGE
#define CSE168_RAW_IMAGE

 enum ImageType {
		RGB,
		RGBA,
		GRAYSCALE,
		HDR
	};

class RawImage {
public:
	RawImage() {}
	RawImage(int w, int h, float* data, ImageType t);
	float* m_rawData;
	int m_width;
	int m_height;
	ImageType m_imageType;

	void loadImage(char* filename);

private:
	void loadHDR(char* filename);
	void loadPPM(char* filename);
    void loadTGA(char* filename);
};

#endif