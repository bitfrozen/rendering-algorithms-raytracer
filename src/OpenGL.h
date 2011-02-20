#ifndef CSE168_OPENGL_H_INCLUDED
#define CSE168_OPENGL_H_INCLUDED

// use the following on Windows or GNU/Linux
#ifndef __APPLE__

#ifndef GLUT_BUILDING_LIB
#define GLUT_BUILDING_LIB
#endif // GLUT_BUILDING_LIB
#ifdef WIN64
#include <freeglut.h>
#else
#include <GL/glut.h>
#endif

#else

// use this on Mac OSX
#include <GLUT/glut.h>

#endif

#endif // CSE168_OPENGL_H_INCLUDED


