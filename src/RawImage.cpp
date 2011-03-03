#include "Image.h"
#include "RawImage.h"
#include "hdrloader.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

RawImage::RawImage(int w, int h, float* data, ImageType t) {
	m_width = w;
	m_height = h;
	m_rawData = data;
	m_imageType = t;
}

void RawImage::loadImage(char* filename) {
	std::string fn(filename); 
	std::string ext = fn.substr(fn.find_last_of(".") + 1);
	if ((ext == "ppm") || (ext == "PPM")) {
		loadPPM(filename);
	} else if ((ext == "hdr") || (ext == "HDR")) {
		loadHDR(filename);
	} else if ((ext == "tga") || (ext == "TGA")) {
		loadTGA(filename);
	}
}


void RawImage::loadHDR(char* filename){
	HDRLoader::load(filename, m_rawData, m_width, m_height);
	m_imageType = HDR;
}
void RawImage::loadPPM(char* filename) {
	const int BUFSIZE = 128;
	FILE* fp;
	unsigned int read;
	unsigned char* rawData;
	char buf[3][BUFSIZE];
	char* retval_fgets;
	size_t retval_sscanf;

	if ( (fp=fopen(filename, "rb")) == NULL)
	{
		std::cerr << "error reading ppm file, could not locate " << filename << std::endl;
		return;
	}

	// Read magic number:
	retval_fgets = fgets(buf[0], BUFSIZE, fp);

	// Read width and height:
	do
	{
		retval_fgets=fgets(buf[0], BUFSIZE, fp);
	} while (buf[0][0] == '#');
	retval_sscanf=sscanf(buf[0], "%s %s", buf[1], buf[2]);
	m_width  = atoi(buf[1]);
	m_height = atoi(buf[2]);

	// Read maxval:
	do
	{
		retval_fgets=fgets(buf[0], BUFSIZE, fp);
	} while (buf[0][0] == '#');

	// Read image data:
	rawData = new unsigned char[m_width * m_height * 3];
	read = fread(rawData, m_width * m_height * 3, 1, fp);
	fclose(fp);
	if (read != 1)
	{
		std::cerr << "error parsing ppm file, incomplete data" << std::endl;
		delete[] rawData;
		m_width = 0;
		m_height = 0;

		return;
	}

	//convert the char array to float for our image.
    m_rawData = new float[m_width * m_height * 3];
	for (int i = m_width * m_height * 3 - 1; i >=0; i--) {
		m_rawData[i] = ((float)rawData[i])/255;
	}
	delete[] rawData;
	m_imageType = RGB;

}
void RawImage::loadTGA(char* filename) {
	FILE *file;
	unsigned char* imageArray;
	int mode,total;

	file = fopen(filename, "rb");
	if (file == NULL) 
	{
		std::cerr << "error parsing tga file, cannot open file" << filename << std::endl;
		return;
	}
	unsigned char type, pixelDepth;
	unsigned char cGarbage;
	short int iGarbage;
	short int width, height;

	// load the header
	fread(&cGarbage, sizeof(unsigned char), 1, file);
	fread(&cGarbage, sizeof(unsigned char), 1, file);
	fread(&type, sizeof(unsigned char), 1, file);
	fread(&iGarbage, sizeof(short int), 1, file);
	fread(&iGarbage, sizeof(short int), 1, file);
	fread(&cGarbage, sizeof(unsigned char), 1, file);
	fread(&iGarbage, sizeof(short int), 1, file);
	fread(&iGarbage, sizeof(short int), 1, file);
	fread(&width, sizeof(short int), 1, file);
	fread(&height, sizeof(short int), 1, file);
	fread(&pixelDepth, sizeof(unsigned char), 1, file);
	fread(&cGarbage, sizeof(unsigned char), 1, file);

	m_width = width;
	m_height = height;
	if ((type != 2) && (type !=3)) 
	{
		std::cerr << "error parsing tga file" << std::endl;
		return;
	}

	// mode equal the number of components for each pixel
	mode = pixelDepth / 8;
	// total is the number of bytes we'll have to read
	total = m_height * m_width * mode;

	// allocate memory for image pixels
	imageArray = new unsigned char[total];

	if (imageArray == NULL) 
	{
		std::cerr << "error parsing tga file, cant make array" << std::endl;
		return;
	}

	fread(imageArray,sizeof(unsigned char),total,file);

	//convert from a char array to a float array
	//Gamma correct loaded data
	m_rawData = new float[total];
	for (int i = total-1; i >= 0; i--) {
		m_rawData[i] = float(Image::gamma_to_linear[imageArray[i]]) / 32768.f;//((float) imageArray[i])/255;
	}
	if (mode == 4)
	{
		for (int i = 3; i < total; i+=4) {
			m_rawData[i] = float(imageArray[i]) / 255.f;//((float) imageArray[i])/255;
		}
	}
	delete[] imageArray;

	if (mode == 1) { //grayscale
		m_imageType = GRAYSCALE;
	} else if (mode == 3) {
		m_imageType = RGB;
	} else if (mode == 4) {
		m_imageType = RGBA;
	}


	// mode=3 or 4 implies that the image is RGB(A). However TGA
	// stores it as BGR(A) so we'll have to swap R and B.
	if (mode >= 3)
	{
		float aux;
		for (int i=0; i < total; i+= mode) 
		{
			aux = m_rawData[i];
			m_rawData[i] = m_rawData[i+2];
			m_rawData[i+2] = aux;
		}
	}
}