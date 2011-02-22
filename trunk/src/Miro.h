#ifndef __MIRO_H__
#define __MIRO_H__

#include <math.h>
#include <stdlib.h>
#include "OpenGL.h"
#include <stdio.h>
#include <iostream>

//#define NO_SSE
#define USE_BINS
#define USE_TRI_PACKETS

#ifdef USE_TRI_PACKETS					// We need SSE for triangle packets...
#ifdef NO_SSE
#undef NO_SSE
#endif
#endif

#define ALIGN_SSE _MM_ALIGN16			// Use this to align variables for SSE use (they need to be 16 byte aligned).
#define ALIGN_64 __declspec(align(64))
typedef unsigned int u_int;

const float MIRO_TMAX = 1e12f;
#ifdef USE_TRI_PACKETS
const int MAX_LEAF_SIZE = 4;
#else
const int MAX_LEAF_SIZE = 2;
#endif
const int NODE_SIZE			= 2;
const int bucket_size		= 32;
const float epsilon			= 0.001f; 
const float PI				= 3.1415926f;
const float piRecip			= 1.0f / PI;
const float DegToRad		= PI/180.0f;
const float RadToDeg		= 180.0f/PI;
const float IntRecip		= (1. / 4294967296.);
const bool use_BVH			= true;
const bool use_Schlick		= false;
const u_int num_samples		= 512;
const u_int NUM_BINS		= 8;
const float INFINITY		= -logf(0);
#ifdef USE_BINS
const bool use_Bins			= true;
#else
const bool use_Bins			= false;
#endif

class Ray;
class HitInfo;

class Object;
class TriangleMesh;
class Instance;

class PointLight;

class Camera;
class Image;
class Scene;
class Material;

extern void ParseFile(FILE* fp);
extern void initOpenGL();
extern Camera* g_camera;
extern Scene* g_scene;
extern Image* g_image;

#endif
