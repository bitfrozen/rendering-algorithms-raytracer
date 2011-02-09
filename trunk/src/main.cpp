#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include "PointLight.h"
#include "Sphere.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"
#include "MiroWindow.h"
#include "Blinn.h"
#include "RawImage.h"
#include "Texture.h"
#include "StoneTexture.h"

void makeMeshObjs(TriangleMesh* mesh, Material* mat);
void makeSpiralScene();
void makeBunnyScene();
void makeBunny20Scene();
void makeEnvironmentMapScene();
void makeStoneFloorScene();
void makeSponzaScene();
void makeTestScene();

int main(int argc, char*argv[])
{
	// create a scene
	makeBunnyScene();
	//makeBunny20Scene();
	//makeStoneFloorScene();
	//makeSponzaScene();
	//makeEnvironmentMapScene();
	//makeTestScene();

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

	// set up the camera
	g_scene->setBGColor(Vector3(0.0,0.0,0.0));
	g_camera->setEye(Vector3(-1.622, 1.361, -1.566));
	g_camera->setLookAt(Vector3(-0.933, 1.57, -0.148));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(45);

	// create and place a point light source
	/*PointLight * light = new PointLight;
	light->setPosition(Vector3(-3, 15, 6));
	light->setColor(Vector3(1, 1, 1));
	light->setWattage(2000);
	g_scene->addLight(light);

	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(-15, 10, -6));
	light2->setColor(Vector3(1, 1, 1));
	light2->setWattage(2000);
	g_scene->addLight(light2);*/

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Mono_Lake.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(0.5f);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(0.5f, 0.5f, 0.5f), Vector3(0.06,0.06,0.06));
	mat->setSpecExp(30.0f);
	mat->setIor(1.58f);
	mat->setReflectAmt(1.0f);
	mat->setRefractAmt(1.0f);
	mat->setEnvMap(hdrTex);
	mat->setEnvExposure(0.5f);
	TriangleMesh *mesh = new TriangleMesh;
	mesh->load("Models/bunny.obj");
	makeMeshObjs(mesh, mat);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeSpiralScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_scene->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addLight(light);

    // create a spiral of spheres
    Material* mat = new Lambert(Vector3(1.0f, 0.0f, 0.0f));
    const int maxI = 200;
    const float a = 0.15f;
    for (int i = 1; i < maxI; ++i)
    {
        float t = i/float(maxI);
        float theta = 4*PI*t;
        float r = a*theta;
        float x = r*cos(theta);
        float y = r*sin(theta);
        float z = 2*(2*PI*a - r);
        Sphere * sphere = new Sphere;
        sphere->setCenter(Vector3(x,y,z));
        sphere->setRadius(r/10);
        sphere->setMaterial(mat);
        g_scene->addObject(sphere);
    }
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void makeBunnyScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(512, 512);

	// set up the camera
	g_scene->setBGColor(Vector3(0.6,0.6,0.85));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(60);

	// create and place a point light source
	PointLight * light = new PointLight;
	light->setPosition(Vector3(-3, 15, 6));
	light->setColor(Vector3(1, 1, 1));
	light->setWattage(2000);
	g_scene->addLight(light);

	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(-15, 10, -6));
	light2->setColor(Vector3(1, 1, 1));
	light2->setWattage(2000);
	g_scene->addLight(light2);

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/Mono_Lake.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.0f);

	// create a spiral of spheres
	Blinn* mat = new Blinn(Vector3(0.5f, 0.5f, 0.5f), Vector3(0.06,0.06,0.06));
	mat->setSpecExp(30.0f);
	mat->setIor(1.56f);
	mat->setReflectAmt(1.0f);
	mat->setRefractAmt(1.0f);
	mat->setEnvMap(hdrTex);
	mat->setEnvExposure(1.0f);
	Blinn* planeMat = new Blinn(Vector3(1.0f, 0.4f, 0.4f), Vector3(0.06,0.06,0.06));
	planeMat->setSpecExp(20.0f);
	planeMat->setSpecAmt(0.1f);
	planeMat->setIor(2.2f);
	planeMat->setReflectAmt(1.0f);
	planeMat->setRefractAmt(0.0f);
	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/bunny.obj");
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
	light2->setWattage(500);
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
	light2->setWattage(500);
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

inline Matrix4x4 translate(float x, float y, float z)
{
	Matrix4x4 m;
	m.setColumn4(Vector4(x, y, z, 1));
	return m;
}


inline Matrix4x4 scale(float x, float y, float z)
{
	Matrix4x4 m;
	m.m11 = x;
	m.m22 = y;
	m.m33 = z;
	return m;
}

// angle is in degrees
inline Matrix4x4 rotate(float angle, float x, float y, float z)
{
	float rad = angle*(PI/180.);

	float x2 = x*x;
	float y2 = y*y;
	float z2 = z*z;
	float c = cos(rad);
	float cinv = 1-c;
	float s = sin(rad);
	float xy = x*y;
	float xz = x*z;
	float yz = y*z;
	float xs = x*s;
	float ys = y*s;
	float zs = z*s;
	float xzcinv = xz*cinv;
	float xycinv = xy*cinv;
	float yzcinv = yz*cinv;

	Matrix4x4 m;
	m.set(x2 + c*(1-x2), xy*cinv+zs, xzcinv - ys, 0,
		xycinv - zs, y2 + c*(1-y2), yzcinv + xs, 0,
		xzcinv + ys, yzcinv - xs, z2 + c*(1-z2), 0,
		0, 0, 0, 1);
	return m;
}

void
makeBunny20Scene()
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
	light->setWattage(1000);
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

	Triangle* t = new Triangle;
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
	light->setWattage(200);
	g_scene->addLight(light);

	Material* material = new Blinn(Vector3(1.0f));
	TriangleMesh * mesh = new TriangleMesh;
	mesh->load("Models/sponza.obj");
	makeMeshObjs(mesh, material);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}

void makeMeshObjs(TriangleMesh* mesh, Material* mat)
{
	int i = mesh->m_numTris-1;
	Triangle* t;
	while (i >= 0)
	{
		t = new Triangle;
		t->setMesh(mesh);
		t->setIndex(i);
		t->setMaterial(mat);
		g_scene->addObject(t);
		i--;
	}
}
