#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include "PointLight.h"
#include "Object.h"
#include "ProxyObject.h"
#include "MBObject.h"
#include "TriangleMesh.h"
#include "Lambert.h"
#include "MiroWindow.h"
#include "Blinn.h"
#include "RawImage.h"
#include "Texture.h"
#include "StoneTexture.h"
#include "Assignment1.h"
#include "Assignment2.h"
#include "BVH.h"

void makeMeshObjs(TriangleMesh* mesh, Material* mat);
void makeMBMeshObjs(TriangleMesh* mesh, TriangleMesh* mesh2, Material* mat);
void makeBunnyScene();
void makeBunnyScene2();
void makeBunny20Scene();
void makeEnvironmentMapScene();
void makeStoneFloorScene();
void makeSponzaScene();
void makeTestScene();
void makeMBTestScene();
void makeProxyTestScene();

int main(int argc, char*argv[])
{
	// create a scene
	//scenes for assignment 1
	//makeBunnyScene();
	//makeTeapotScene();
	//makeSphereScene();

	//scenes that show other functionality
	//makeBunnyScene2();
	//makeBunny20Scene();
	//makeStoneFloorScene();
	//makeSponzaScene();
	//makeEnvironmentMapScene();
	//makeTestScene();
	//makeMBTestScene();
	//makeProxyTestScene();

	//assignment 2
	makePathTracingScene3();
	//makeTeapotScene2();
	//makeBunny1Scene2();
	//makeBunny20Scene2();
	//makeSponzaScene2();

	MiroWindow miro(&argc, argv);
	miro.mainLoop();

	return 0; // never executed
}

void makeTestScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	g_scene->m_pathTrace = true;
	g_scene->m_numPaths = 1;
	g_scene->m_maxBounces = 5;
	g_scene->m_minSubdivs = 1;
	g_scene->m_maxSubdivs = 6;
	g_scene->setNoise(0.03f);

	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(60);
	g_camera->m_aperture = 0.1f;
	g_camera->m_focusPlane = 3.0f;

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Topanga_Forest_B_3k.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(0.09, 0.094, 0.1));
	mat->setSpecExp(30.0f);
	mat->setSpecAmt(0);
	mat->setIor(6.0f);
	mat->setReflectAmt(0.90f);
	mat->setRefractAmt(0.0f);
	mat->setReflectGloss(0.98f);

	// create a spiral of spheres
	Blinn* mat2 = new Blinn(Vector3(0.4, 0.05, 0.02), Vector3(0), Vector3(1.0, 0.5, 0.0));
	mat2->setSpecExp(30.0f);
	mat2->setSpecAmt(0);
	mat2->setIor(8.0f);
	mat2->setReflectAmt(0.90f);
	mat2->setRefractAmt(0.0f);
	mat2->setReflectGloss(0.98f);

	Blinn* planeMat = new Blinn(Vector3(1.0));
	planeMat->setSpecExp(20.0f);
	planeMat->setSpecAmt(0);
	planeMat->setIor(1.6f);
	planeMat->setReflectAmt(1.0f);
	planeMat->setRefractAmt(0.0f);
	planeMat->setReflectGloss(0.97f);

	Matrix4x4 xform;
	xform *= scale(0.8, 0.8, 0.8);
	xform *= rotate(30, 0, 1, 0);
	xform *= translate(-5, 0, 0);

	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *mesh2 = new TriangleMesh;
	TriangleMesh *plane = new TriangleMesh;
	mesh->load("Models/bunny.obj");
	mesh2->load("Models/bunny.obj", xform);
	plane->load("Models/plane.obj");
	makeMeshObjs(mesh, mat);
	makeMeshObjs(mesh2, mat2);
	makeMeshObjs(plane, planeMat);

	Object* t = new Object;

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeMBTestScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	g_scene->m_pathTrace = false;
	g_scene->m_numPaths = 1;
	g_scene->m_maxBounces = 5;
	g_scene->m_minSubdivs = 3;
	g_scene->m_maxSubdivs = 16;
	g_scene->setNoise(0.002f);

	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-5.536, 9.552, 15.815));
	g_camera->setViewDir(Vector3(-0.034, -0.009, -0.999));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(75);
	g_camera->m_aperture = 0.001f;
	g_camera->m_focusPlane = 3.0f;
	g_camera->setShutterSpeed(0.2f);

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Topanga_Forest_B_3k.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);

	// create and place a point light source
	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(5, 5, 5));
	light2->setColor(Vector3(1, 1, 1));
	light2->setPower(500);
	g_scene->addLight(light2);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(0.09, 0.094, 0.1));
	mat->setSpecExp(30.0f);
	mat->setSpecAmt(0);
	mat->setIor(6.0f);
	mat->setReflectAmt(0.90f);
	mat->setRefractAmt(0.0f);
	mat->setReflectGloss(0.98f);

	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *mesh2 = new TriangleMesh;
	mesh->load("Models/bulletMB_01.obj");
	mesh2->load("Models/bulletMB_02.obj");
	makeMBMeshObjs(mesh, mesh2, mat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeProxyGrid(Objects* o, BVH* b);

void makeProxyTestScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	g_scene->m_pathTrace = true;
	g_scene->m_numPaths = 2;
	g_scene->m_maxBounces = 3;
	g_scene->m_minSubdivs = 3;
	g_scene->m_maxSubdivs = 4;
	g_scene->setNoise(0.01f);

	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(60);

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Topanga_Forest_B_3k.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(1));//0.09, 0.094, 0.1));
	mat->setSpecExp(30.0f);
	mat->setSpecAmt(0);
	mat->setIor(6.0f);
	mat->setReflectAmt(0.0f);
	mat->setRefractAmt(0.0f);
	mat->setReflectGloss(0.98f);

	Blinn* planeMat = new Blinn(Vector3(1.0));
	planeMat->setSpecExp(20.0f);
	planeMat->setSpecAmt(0);
	planeMat->setIor(1.6f);
	planeMat->setReflectAmt(0.0f);
	planeMat->setRefractAmt(0.0f);
	planeMat->setReflectGloss(1.f);

	TriangleMesh *mesh = new TriangleMesh;
	//TriangleMesh *plane = new TriangleMesh;

	mesh->load("Models/bunny.obj");
	//plane->load("Models/plane.obj");

	Matrix4x4 m = Matrix4x4();
	m.translate(0,0,2);

	BVH* b = new BVH;
	Objects* o = new Objects;
	ProxyObject::setupProxy(mesh, mat, o, b);

	makeProxyGrid(o, b);

	//makeMeshObjs(plane, planeMat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeProxyGrid(Objects* o, BVH* b)
{
	for (int i = 0; i <= 1000; i++)
	{
		for (int j = 0; j <= 1000; j++)
		{
			Matrix4x4 m = Matrix4x4();
			m.translate(2*(i-500), 0, 2*(j-500));
			ProxyObject* po = new ProxyObject(o, b, m);
			g_scene->addObject(po);
		}
	}
}

void makeBunnyScene2()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	g_scene->m_pathTrace = false;
	g_scene->m_numPaths = 1;
	g_scene->m_maxBounces = 4;
	g_scene->setNoise(0.01f);
	g_scene->setMinSubdivs(1);
	g_scene->setMaxSubdivs(6);

	// set up the camera
	g_scene->setBGColor(Vector3(1.0,0.0,0.0));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0.65, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(45);

	// create and place a point light source
	PointLight * light = new PointLight;
	light->setPosition(Vector3(-3, 15, 6));
	light->setColor(Vector3(1, 1, 1));
	light->setPower(2000);
	g_scene->addLight(light);

	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(-15, 10, -6));
	light2->setColor(Vector3(1, 1, 1));
	light2->setPower(2000);
	g_scene->addLight(light2);

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Mono_Lake.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(0.7f, 0.1f, 0.05f), Vector3(0));
	mat->setSpecExp(30.0f);
	mat->setIor(2.2f);
	mat->setReflectAmt(1.0f);
	mat->setRefractAmt(1.0f);
	mat->setEnvMap(hdrTex);
	mat->setEnvExposure(1.0f);

	Blinn* planeMat = new Blinn(Vector3(1.0f, 0.4f, 0.4f), Vector3(0));
	planeMat->setSpecExp(20.0f);
	planeMat->setSpecAmt(0.1f);
	planeMat->setIor(2.2f);
	planeMat->setReflectAmt(1.0f);
	planeMat->setRefractAmt(0.0f);
	planeMat->setReflectGloss(0.95f);

	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/sphere2.obj");
	planeMesh->load("Models/plane.obj");
	makeMeshObjs(mesh, mat);
	makeMeshObjs(planeMesh, planeMat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeStoneFloorScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	// set up the camera
	g_scene->setBGColor(Vector3(0.6,0.6,0.85));
	g_camera->setEye(Vector3(0, 5, 10));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(90);

	// create and place a point light source
	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(5, 5, 5));
	light2->setColor(Vector3(1, 1, 1));
	light2->setPower(500);
	g_scene->addLight(light2);

	//some texture stuff
	Material* stoneMat = new Lambert(Vector3(0.5f, 0.5f, 0.5f), Vector3(0.06,0.06,0.06));
	StoneTexture* stoneTex = new StoneTexture(100);
	stoneMat->m_texture = stoneTex;
	//end

	// create sphere
	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/sphere2.obj");
	planeMesh->load("Models/plane.obj");
	makeMeshObjs(mesh, stoneMat);
	makeMeshObjs(planeMesh, stoneMat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeEnvironmentMapScene() {
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	// set up the camera
	g_scene->setBGColor(Vector3(0.6,0.6,0.85));
	g_camera->setEye(Vector3(0, 5, -10));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(90);

	// create and place a point light source
	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(-5, 10, -5));
	light2->setColor(Vector3(1, 1, 1));
	light2->setPower(500);
	g_scene->addLight(light2);

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Mono_Lake.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(2.0f);

	//generate stone texture
	Material* stoneMat = new Lambert(Vector3(0.5f, 0.5f, 0.5f), Vector3(0.06,0.06,0.06));
	StoneTexture* stoneTex = new StoneTexture(150);
	stoneMat->m_texture = stoneTex;
	//end

	// create sphere
	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/sphere2.obj");
	planeMesh->load("Models/plane.obj");
	makeMeshObjs(mesh, stoneMat);
	makeMeshObjs(planeMesh, stoneMat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeBunny20Scene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	// set up the camera
	g_scene->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
	g_camera->setEye(Vector3(0, 5, 15));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(45);

	// create and place a point light source
	PointLight * light = new PointLight;
	light->setPosition(Vector3(10, 20, 10));
	light->setColor(Vector3(1, 1, 1));
	light->setPower(1000);
	g_scene->addLight(light);

	TriangleMesh * mesh;
	Blinn* material = new Blinn(Vector3(1.0f));

	Blinn* mat2 = new Blinn(Vector3(1.0f));
	Matrix4x4 xform;
	Matrix4x4 xform2;
	xform2 *= rotate(110, 0, 1, 0);
	xform2 *= scale(.6, 1, 1.1);


	// bunny 1
	xform.setIdentity();
	xform *= scale(0.3, 2.0, 0.7);
	xform *= translate(-1, .4, .3);
	xform *= rotate(25, .3, .1, .6);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 2
	xform.setIdentity();
	xform *= scale(.6, 1.2, .9);
	xform *= translate(7.6, .8, .6);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 3
	xform.setIdentity();
	xform *= translate(.7, 0, -2);
	xform *= rotate(120, 0, .6, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 4
	xform.setIdentity();
	xform *= translate(3.6, 3, -1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 5
	xform.setIdentity();
	xform *= translate(-2.4, 2, 3);
	xform *= scale(1, .8, 2);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 6
	xform.setIdentity();
	xform *= translate(5.5, -.5, 1);
	xform *= scale(1, 2, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 7
	xform.setIdentity();
	xform *= rotate(15, 0, 0, 1);
	xform *= translate(-4, -.5, -6);
	xform *= scale(1, 2, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 8
	xform.setIdentity();
	xform *= rotate(60, 0, 1, 0);
	xform *= translate(5, .1, 3);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 9
	xform.setIdentity();
	xform *= translate(-3, .4, 6);
	xform *= rotate(-30, 0, 1, 0);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 10
	xform.setIdentity();
	xform *= translate(3, 0.5, -2);
	xform *= rotate(180, 0, 1, 0);
	xform *= scale(1.5, 1.5, 1.5);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 11
	xform = xform2;
	xform *= scale(0.3, 2.0, 0.7);
	xform *= translate(-1, .4, .3);
	xform *= rotate(25, .3, .1, .6);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 12
	xform = xform2;
	xform *= scale(.6, 1.2, .9);
	xform *= translate(7.6, .8, .6);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 13
	xform = xform2;
	xform *= translate(.7, 0, -2);
	xform *= rotate(120, 0, .6, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 14
	xform = xform2;
	xform *= translate(3.6, 3, -1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 15
	xform = xform2;
	xform *= translate(-2.4, 2, 3);
	xform *= scale(1, .8, 2);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 16
	xform = xform2;
	xform *= translate(5.5, -.5, 1);
	xform *= scale(1, 2, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 17
	xform = xform2;
	xform *= rotate(15, 0, 0, 1);
	xform *= translate(-4, -.5, -6);
	xform *= scale(1, 2, 1);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 18
	xform = xform2;
	xform *= rotate(60, 0, 1, 0);
	xform *= translate(5, .1, 3);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 19
	xform = xform2;
	xform *= translate(-3, .4, 6);
	xform *= rotate(-30, 0, 1, 0);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// bunny 20
	xform = xform2;
	xform *= translate(3, 0.5, -2);
	xform *= rotate(180, 0, 1, 0);
	xform *= scale(1.5, 1.5, 1.5);
	mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj", xform);
	makeMeshObjs(mesh, material);

	// create the floor triangle
	mesh = new TriangleMesh;
	mesh->createSingleTriangle();
	mesh->setV1(Vector3(-100, 0, -100));
	mesh->setV2(Vector3(   0, 0,  100));
	mesh->setV3(Vector3( 100, 0, -100));
	mesh->setN1(Vector3(0, 1, 0));
	mesh->setN2(Vector3(0, 1, 0));
	mesh->setN3(Vector3(0, 1, 0));

	Object* t = new Object;
	t->setIndex(0);
	t->setMesh(mesh);
	t->setMaterial(mat2); 
	g_scene->addObject(t);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void
makeSponzaScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	// set up the camera
	g_scene->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
	g_camera->setEye(Vector3(8, 1.5, 1));
	g_camera->setLookAt(Vector3(0, 2.5, -1));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(55);

	// create and place a point light source
	PointLight * light = new PointLight;
	light->setPosition(Vector3(0, 10.0, 0));
	light->setColor(Vector3(1, 1, 1));
	light->setPower(200);
	g_scene->addLight(light);

	Material* material = new Blinn(Vector3(1.0f));
	TriangleMesh* mesh = new TriangleMesh;
	mesh->load("Models/sponza.obj");
	makeMeshObjs(mesh, material);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

// Insert all objects into the scene in a single memory block
void makeMeshObjs(TriangleMesh* mesh, Material* mat)
{
	int numObjs = mesh->m_numTris;

	Object* t = new Object[numObjs];//(Object*)_aligned_malloc(sizeof(Object)*numObjs, 16);
	
	int i = numObjs-1;
	while (i >= 0)
	{
		t[i].setMesh(mesh);
		t[i].setIndex(i);
		t[i].setMaterial(mat);
		g_scene->addObject(&t[i]);
		i--;
	}
}

// Insert all MotionBlurred objects into the scene in a single memory block
void makeMBMeshObjs(TriangleMesh* mesh, TriangleMesh* mesh2, Material* mat)
{
	int numObjs = mesh->m_numTris;

	MBObject* t = new MBObject[numObjs];//(MBObject*)_aligned_malloc(sizeof(MBObject)*numObjs, 16);

	int i = numObjs-1;
	while (i >= 0)
	{
		t[i].setMesh(mesh);
		t[i].setMeshT2(mesh2);
		t[i].setIndex(i);
		t[i].setMaterial(mat);
		g_scene->addObject(&t[i]);
		i--;
	}
}