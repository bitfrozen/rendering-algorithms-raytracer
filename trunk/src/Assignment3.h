#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"

#include "PointLight.h"
#include "TriangleMesh.h"
#include "Lambert.h"
#include "Material.h"
#include "Blinn.h"
#include "RectangleLight.h"
#include "ProxyObject.h"

void makeAlphaTest();

void testSceneTree();

void
makeAlphaTest()
{
    g_camera = (Camera*)_aligned_malloc(sizeof(Camera), 16);
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(256, 256);

	g_scene->m_pathTrace = true;
	g_scene->m_numPaths = 20;
	g_scene->m_maxBounces = 20;
	g_scene->m_minSubdivs = 1;
	g_scene->m_maxSubdivs = 4;
	g_scene->setNoise(0.01f);
    
    // set up the camera
    g_scene->setBGColor(Vector3(0.0f));
    g_camera->setEye(Vector3(0, 3, 6));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);
	g_camera->m_aperture = 0.001f;
	g_camera->m_focusPlane = 4.0f;

    // create and place a point light source
    /*PointLight * light = new PointLight;
    light->setPosition(Vector3(0, 10, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setPower(700);
    g_scene->addLight(light);
*/
	PointLight * light2 = new PointLight;
    light2->setPosition(Vector3(-10, -10, -10));
    light2->setColor(Vector3(1, 1, 1));
    light2->setPower(4000);
    g_scene->addLight(light2);

    Material* leafMaterial = new Blinn(Vector3(1.0f));
	leafMaterial->setTranslucency(.0f);
	RawImage* leafRI = new RawImage();
	leafRI->loadImage("Textures/Tree_03_Leaves.tga");
	Texture* leafTex = new Texture(leafRI);
	leafMaterial->setColorMap(leafTex);
	leafMaterial->setAlphaMap(leafTex);

//     TriangleMesh * leafMesh = new TriangleMesh;
//     leafMesh->load("Models/leaf_test.obj");
//     makeMeshObjs(leafMesh, leafMaterial);

	Material* leafMaterial2 = new Blinn(Vector3(1.0f));
	leafMaterial2->setTranslucency(.90f);
	leafMaterial2->setColorMap(leafTex);
	leafMaterial2->setAlphaMap(leafTex);
	Matrix4x4 xform;
	xform = translate(-2, 0, 0);
	TriangleMesh* mesh = new TriangleMesh;
    mesh->load("Models/leaf_test.obj", xform);
    makeMeshObjs(mesh, leafMaterial2);

	xform = translate(-1,.5,0);
	TriangleMesh* leaf3 = new TriangleMesh;
	leaf3->load("Models/leaf_test.obj", xform);
	makeMeshObjs(leaf3, leafMaterial2);
    
    // create the floor triangle
	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Topanga_Forest_B_3k.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);
	g_scene->setSampleEnv(true);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}