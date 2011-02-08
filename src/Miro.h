#ifndef __MIRO_H__
#define __MIRO_H__

//#define NO_SSE
#define ALIGN_SSE __declspec(align(16))			// Use this to align variables for SSE use (they need to be 16 byte aligned).
typedef unsigned int u_int;

// #ifndef max
// #define max(a,b) ((a>b)?a:b)
// #endif

const float MIRO_TMAX = 1e12f;
const int MAX_LEAF_SIZE = 8;
const int NODE_SIZE = 2;
const float epsilon   = 0.001f; 
const float PI = 3.1415926f;
const float DegToRad = PI/180.0f;
const float RadToDeg = 180.0f/PI;
const float Ttri = 1.0f;
const float Tbox = 1.0f;
const bool use_BVH = true;

#include <stdlib.h>
#include "OpenGL.h"
#include <stdio.h>
#include <iostream>

class Ray;
class HitInfo;

class Object;
class Sphere;
class Triangle;
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
