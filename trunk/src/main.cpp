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
#include "TextureGenerator.h"
#include "Texture.h"

void
makeSpiralScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
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

void
makeBunnyScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(256, 256);

	// set up the camera
	g_camera->setBGColor(Vector3(0.6,0.6,0.85));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(45);

	// create and place a point light source
	PointLight * light = new PointLight;
	light->setPosition(Vector3(-3, 15, 6));
	light->setColor(Vector3(1, 1, 1));
	light->setWattage(2000);
//	g_scene->addLight(light);

	PointLight * light2 = new PointLight;
	//light2->setPosition(Vector3(-15, 10, -6));
	light2->setPosition(Vector3(-5, 2, 3));
	light2->setColor(Vector3(1, 1, 1));
	light2->setWattage(500);
	g_scene->addLight(light2);

	// create a spiral of spheres
	Material* mat = new Blinn(Vector3(0.1f, 0.5f, 0.2f), Vector3(0.06,0.06,0.06));
	//((Blinn *)mat)->setSpecExp(10.0f);
	//((Blinn *)mat)->setReflectAmt(0.4f);
	//((Blinn *)mat)->setRefractAmt(0.6f);
	Material* planeMat = new Blinn(Vector3(0.4f, 0.4f, 0.4f), Vector3(0.06,0.06,0.06));
	((Blinn *)planeMat)->setSpecExp(20.0f);
	((Blinn *)planeMat)->setSpecAmt(0.1f);
	//((Blinn *)planeMat)->setReflectAmt(0.4f);
	//((Blinn *)planeMat)->setRefractAmt(0.6f);
	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/plane.obj");
//	planeMesh->load("Models/plane.obj");
	Triangle *bunny = new Triangle;
	Triangle *plane = new Triangle;
	bunny->setMesh(mesh);
	bunny->setMaterial(planeMat);
	plane->setMesh(planeMesh);
	plane->setMaterial(planeMat);
	g_scene->addObject(bunny);
//	g_scene->addObject(plane);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}


void
makeStoneFloorScene()
{
	g_camera = new Camera;
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(256, 256);

	// set up the camera
	g_camera->setBGColor(Vector3(0.6,0.6,0.85));
	g_camera->setEye(Vector3(0, 5, 10));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(45);

	// create and place a point light source
	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(5, 5, 5));
	//light2->setPosition(Vector3(-5, 2, 3));
	light2->setColor(Vector3(1, 1, 1));
	light2->setWattage(500);
	g_scene->addLight(light2);

	//some texture stuff
	Material* cellMat = new Lambert(Vector3(0.5f, 0.5f, 0.5f), Vector3(0.06,0.06,0.06));
	TextureGenerator* tg = new TextureGenerator();
	Texture* cellTex = tg->generateStoneTexture(300, 300, 200);
	delete tg;
	cellMat->m_texture = cellTex;
	//end

	// create sphere
	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *planeMesh = new TriangleMesh;
	mesh->load("Models/sphere2.obj");
	planeMesh->load("Models/plane.obj");
	Triangle *bunny = new Triangle;
	Triangle *plane = new Triangle;
	bunny->setMesh(mesh);
	bunny->setMaterial(cellMat);
	plane->setMesh(planeMesh);
	plane->setMaterial(cellMat);
	g_scene->addObject(bunny);
	g_scene->addObject(plane);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}
int
main(int argc, char*argv[])
{
	// create a scene
	makeBunnyScene()
    //makeStoneFloorScene();

    MiroWindow miro(&argc, argv);
    miro.mainLoop();

    return 0; // never executed
}

