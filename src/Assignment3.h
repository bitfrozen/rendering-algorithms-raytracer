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

void testDispersion();



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

void
testDispersion()
{
    g_camera = (Camera*)_aligned_malloc(sizeof(Camera), 16);
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(256, 256);

	g_scene->m_pathTrace = true;
	g_scene->m_numPaths = 20;
	g_scene->m_maxBounces = 6;
	g_scene->m_minSubdivs = 1;
	g_scene->m_maxSubdivs = 4;
	g_scene->setNoise(0.01f);
    
    // set up the camera
    g_scene->setBGColor(Vector3(0.0f));
    g_camera->setEye(Vector3(0, 3, 6));
    g_camera->setLookAt(Vector3(0, 2, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);
	g_camera->m_aperture = 0.001f;
	g_camera->m_focusPlane = 4.0f;

		//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/sky.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.f);
	g_scene->setSampleEnv(false);
	
	//make a raw image from hdr file
	RawImage* hdrImage2 = new RawImage();
	hdrImage2->loadImage("Images/sky.hdr");
	Texture* hdrTex2 = new Texture(hdrImage2);

	DomeLight* domeLight = new DomeLight; 
	domeLight->setTexture(hdrTex2);
	domeLight->setPower(0.15f);
	domeLight->setSamples(6);
	g_scene->addLight(domeLight);

	//make a raw image from hdr file
	RawImage* hdrImage3 = new RawImage();
	hdrImage3->loadImage("Images/Topanga_Forest_B_3k.hdr");
	Texture* hdrTex3 = new Texture(hdrImage3);
	g_scene->setEnvMap(hdrTex3);
	g_scene->setEnvExposure(1.0f);

    // create and place a point light source
    /*PointLight * light = new PointLight;
    light->setPosition(Vector3(0, 10, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setPower(700);
    g_scene->addLight(light);
*/
	/*PointLight * light2 = new PointLight;
    light2->setPosition(Vector3(-10, -10, -10));
    light2->setColor(Vector3(1, 1, 1));
    light2->setPower(.04);
    g_scene->addLight(light2);
*/

	Blinn* mat = new Blinn(Vector3(0.0f, 0.5f, 0.5f), Vector3(0.0,0.0,0.0));
	mat->setSpecExp(30.0f);
	mat->setIor(1.608f);
	mat->setReflectAmt(1.0f);
	mat->setRefractAmt(1.0f);
	mat->m_disperse = false;

	Blinn* mat2 = new Blinn(Vector3(0.0f, 0.5f, 0.5f), Vector3(0.0,0.0,0.0));
	mat2->setSpecExp(30.0f);
	mat2->setIor(1.608f);
	mat2->setReflectAmt(1.0f);
	mat2->setRefractAmt(1.0f);
	mat2->m_disperse = true;

    
	Matrix4x4 xform;
	xform = translate(0, 0, 0);
	TriangleMesh* mesh = new TriangleMesh;
    mesh->load("Models/sphere2.obj", xform);
    makeMeshObjs(mesh, mat2);

	/*xform = translate(-4,0,0);
	TriangleMesh* mesh2 = new TriangleMesh;
	mesh2->load("Models/sphere2.obj", xform);
	makeMeshObjs(mesh2, mat2);
    */
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void testSceneTree() {

	g_camera = (Camera*)_aligned_malloc(sizeof(Camera), 16);
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(256, 256);

	g_scene->m_pathTrace = true;
	g_scene->m_numPaths = 2;
	g_scene->m_maxBounces = 3;
	g_scene->m_minSubdivs = 1;
	g_scene->m_maxSubdivs = 4;
	g_scene->setNoise(0.01f);

	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-5, 4, 3));
	g_camera->setLookAt(Vector3(0, 0, 0));
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(60);
	g_camera->m_aperture = 0.001f;
	g_camera->m_focusPlane = 4.0f;

	//make a raw image from hdr file
	RawImage* hdrImage = new RawImage();
	hdrImage->loadImage("Images/sky.hdr");
	Texture* hdrTex = new Texture(hdrImage);
	g_scene->setEnvMap(hdrTex);
	g_scene->setEnvExposure(1.f);
	g_scene->setSampleEnv(false);
	
	//make a raw image from hdr file
	RawImage* hdrImage2 = new RawImage();
	hdrImage2->loadImage("Images/sky.hdr");
	Texture* hdrTex2 = new Texture(hdrImage2);

	DomeLight* domeLight = new DomeLight; 
	domeLight->setTexture(hdrTex2);
	domeLight->setPower(0.15f);
	domeLight->setSamples(6);
	g_scene->addLight(domeLight);

	/*// create and place a point light source
	PointLight * light2 = new PointLight;
	light2->setPosition(Vector3(-5, 5, 0));
	light2->setColor(Vector3(1, 1, 1));
	light2->setPower(500);
	g_scene->addLight(light2);*/

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

	Blinn* planeMat = new Blinn(Vector3(0.2f));
	planeMat->setSpecExp(10.0f);
	planeMat->setSpecAmt(0);
	planeMat->setIor(10.f);
	planeMat->setReflectAmt(0.0f);
	planeMat->setRefractAmt(0.0f);
	planeMat->setReflectGloss(1.0f);

	RawImage* barkImg = new RawImage();
	barkImg->loadImage("Textures/bark_COLOR.tga");
	Texture* barkTex = new Texture(barkImg);

	RawImage* barkNrmImg = new RawImage();
	barkNrmImg->loadImage("Textures/bark_NRM.tga");
	Texture* barkNrmTex = new Texture(barkNrmImg);

	RawImage* barkSpecImg = new RawImage();
	barkSpecImg->loadImage("Textures/bark_SPEC.tga");
	Texture* barkSpecTex = new Texture(barkSpecImg);

	Blinn* barkMtl = new Blinn(Vector3(1.0));
	barkMtl->setSpecExp(10.0f);
	barkMtl->setSpecAmt(0.5);
	barkMtl->setIor(1.6f);
	barkMtl->setReflectAmt(0.0f);
	barkMtl->setRefractAmt(0.0f);
	barkMtl->setReflectGloss(1.0f);
	barkMtl->setColorMap(barkTex);
	barkMtl->setNormalMap(barkNrmTex);
	barkMtl->setSpecularMap(barkSpecTex);

	RawImage* leavesImg = new RawImage();
	leavesImg->loadImage("Textures/testTreeLeaves.tga");
	Texture* leavesTex = new Texture(leavesImg);

	Blinn* leavesMtl = new Blinn(Vector3(0.5));
	leavesMtl->setSpecExp(10.0f);
	leavesMtl->setSpecAmt(0.5);
	leavesMtl->setColorMap(leavesTex);
	leavesMtl->setAlphaMap(leavesTex);
	leavesMtl->setTranslucency(.9f);

	RawImage* grassImg = new RawImage();
	grassImg->loadImage("Textures/grassblade2.tga");
	Texture* grassTex = new Texture(grassImg);

	Blinn* grassMtl = new Blinn(Vector3(0.5));
	grassMtl->setSpecExp(20.0f);
	grassMtl->setSpecAmt(0.8);
	grassMtl->setColorMap(grassTex);

	Matrix4x4 xform;
	xform.rotate(30, 0, 1, 0);
	xform.translate(-3.5, 0.05, -2);
	xform.scale(0.8, 0.8, 0.8);
	
// 	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *mesh2 = new TriangleMesh;
	TriangleMesh* tree = new TriangleMesh;
	TriangleMesh* leaves = new TriangleMesh;
	TriangleMesh *plane = new TriangleMesh;
	TriangleMesh *grass = new TriangleMesh;
//	mesh->load("Models/bunny.obj");
	mesh2->load("Models/bunny.obj", xform);
	tree->load("Models/testTree.obj");
	leaves->load("Models/testTreeLeaves.obj");
	plane->load("Models/plane.obj");
	grass->load("Models/testGrass2.obj");
//  	makeMeshObjs(mesh, mat);
	makeMeshObjs(mesh2, mat);
	makeMeshObjs(tree, barkMtl);
	makeMeshObjs(leaves, leavesMtl);
	makeMeshObjs(plane, planeMat);
 
 	BVH* b = new BVH;
 	Objects* o = new Objects;
 	//ProxyObject::setupProxy(grass, grassMtl, o, b);

	//makeProxyGrid(o, b);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}