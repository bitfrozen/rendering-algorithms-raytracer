#include <stdio.h>
#include <stdlib.h>
#include "Miro.h"
#include "Camera.h"
#include "Image.h"
#include "Scene.h"
#include "Console.h" 
#include "OpenGL.h"
#include <time.h>

Camera * g_camera = 0;

const float HalfDegToRad = DegToRad/2.0f;

Camera::Camera() :
    m_renderer(RENDER_OPENGL),
    m_eye(0,0,0),
    m_viewDir(0,0,-1),
    m_up(0,1,0),
    m_lookAt(FLT_MAX, FLT_MAX, FLT_MAX),
    m_fov((45.)*(PI/180.))
{
    calcLookAt();
}

Camera::~Camera()
{
}

void Camera::click(Scene* pScene, Image* pImage)
{
    calcLookAt();
    static bool firstRayTrace = false;

    if (m_renderer == RENDER_OPENGL)
    {
        glDrawBuffer(GL_BACK);
        pScene->openGL(this);
        firstRayTrace = true;
    }
    else if (m_renderer == RENDER_RAYTRACE)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDrawBuffer(GL_FRONT);
        if (firstRayTrace)
        {
            pImage->clear(Vector3(0));
            pScene->raytraceImage(this, g_image);
            firstRayTrace = false;
        }        
        g_image->draw();
    }
}

void Camera::calcLookAt()
{
    // this is true when a "lookat" is not used in the config file
    if (m_lookAt.x != FLT_MAX)
    {
        setLookAt(m_lookAt);
        m_lookAt.set(FLT_MAX, FLT_MAX, FLT_MAX);
    }
}

void Camera::drawGL()
{
    // set up the screen with our camera parameters
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov(), g_image->width()/(float)g_image->height(),
                   0.01, 10000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vector3 vCenter = eye() + viewDir();
    gluLookAt(eye().x, eye().y, eye().z,
              vCenter.x, vCenter.y, vCenter.z,
              up().x, up().y, up().z);
}

Ray Camera::eyeRay(int x, int y, int imageWidth, int imageHeight)
{
    // first compute the camera coordinate system 
    // ------------------------------------------
    // wDir = e - (e+m_viewDir) = -m_vView
    const Vector3 wDir = Vector3(-m_viewDir).normalize(); 
    const Vector3 uDir = cross(m_up, wDir).normalize(); 
    const Vector3 vDir = cross(wDir, uDir);    

    // next find the corners of the image plane in camera space
    // --------------------------------------------------------
    const float aspectRatio = (float)imageWidth/(float)imageHeight; 

    const float top     = tan(m_fov*HalfDegToRad); 
    const float right   = aspectRatio*top; 

    const float bottom  = -top; 
    const float left    = -right; 

    // transform x and y into camera space 
    // -----------------------------------

    const float imPlaneUPos = left   + (right - left)*(((float)x+0.5f)/(float)imageWidth); 
    const float imPlaneVPos = bottom + (top - bottom)*(((float)y+0.5f)/(float)imageHeight);

    return Ray(m_eye, (imPlaneUPos*uDir + imPlaneVPos*vDir - wDir).normalize());
}

Ray Camera::eyeRayRandom(int x, int y, int imageWidth, int imageHeight)
{
    // first compute the camera coordinate system 
    // ------------------------------------------

    // wDir = e - (e+m_viewDir) = -m_vView
    const Vector3 wDir = Vector3(-m_viewDir).normalize(); 
    const Vector3 uDir = cross(m_up, wDir).normalize(); 
    const Vector3 vDir = cross(wDir, uDir);    

    // next find the corners of the image plane in camera space
    // --------------------------------------------------------

    const float aspectRatio = (float)imageWidth/(float)imageHeight; 

    const float top     = tan(m_fov*HalfDegToRad); 
    const float right   = aspectRatio*top; 

    const float bottom  = -top; 
    const float left    = -right; 

    // transform x and y into camera space 
    // -----------------------------------
    //add a randomness so that the ray will be sent somewhere inside the pixel,
	// and not necessarily through the center of the pixel
	if (Scene::randsIdx >= 990000)
	{
		Scene::randsIdx = 0;
		Scene::genRands();
	}
	float urand = Scene::rands[Scene::randsIdx++];
	float vrand = Scene::rands[Scene::randsIdx++];

    const float imPlaneUPos = left   + (right - left)*(((float)x+urand)/(float)imageWidth); 
    const float imPlaneVPos = bottom + (top - bottom)*(((float)y+vrand)/(float)imageHeight); 

    return Ray(m_eye, (imPlaneUPos*uDir + imPlaneVPos*vDir - wDir).normalize());
}

Ray Camera::eyeRayRandomDOF(int x, int y, int imageWidth, int imageHeight)
{
	// first compute the camera coordinate system 
	// ------------------------------------------

	// wDir = e - (e+m_viewDir) = -m_vView
	const Vector3 wDir = Vector3(-m_viewDir).normalize(); 
	const Vector3 uDir = cross(m_up, wDir).normalize(); 
	const Vector3 vDir = cross(wDir, uDir);

	// next find the corners of the image plane in camera space
	// --------------------------------------------------------

	const float aspectRatio = (float)imageWidth/(float)imageHeight; 

	const float top     = tan(m_fov*HalfDegToRad); 
	const float right   = aspectRatio*top; 

	const float bottom  = -top; 
	const float left    = -right; 

	// transform x and y into camera space 
	// -----------------------------------
	//add a randomness so that the ray will be sent somewhere inside the pixel,
	// and not necessarily through the center of the pixel
	float urand = Scene::getRand();
	float vrand = Scene::getRand();

	const float imPlaneUPos = left   + (right - left)*(((float)x+urand)/(float)imageWidth);
	const float imPlaneVPos = bottom + (top - bottom)*(((float)y+vrand)/(float)imageHeight);

	const Vector3 focalPoint = (imPlaneUPos*uDir + imPlaneVPos*vDir - wDir).normalized() * m_focusPlane + m_eye;

	// Rejection sample a disc
	do 
	{
		urand = 1.0 - 2*Scene::getRand();
		vrand = 1.0 - 2*Scene::getRand();
	} while (urand*urand + vrand*vrand > 1.0f);

	const Vector3 origin = m_aperture*(urand*uDir + vrand*vDir) + m_eye;
	const Vector3 direction = (focalPoint - origin).normalized();

	return Ray(origin, direction);
}