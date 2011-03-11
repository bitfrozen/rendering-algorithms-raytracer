#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include "PointLight.h"
#include "DomeLight.h"
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
#include "BVH.h"

void makeFinalScene();
void makeMeshObjs(TriangleMesh* mesh, Material* mat);
void makeMBMeshObjs(TriangleMesh* mesh, TriangleMesh* mesh2, Material* mat);
void makeProxyGrid(Objects* o, BVH* b, int n);

int main(int argc, char*argv[])
{
	// create a scene
	makeFinalScene();

	MiroWindow miro(&argc, argv);
	miro.mainLoop();

	return 0; // never executed
}

void makeProxyGrid(Objects* o, BVH* b, int n)
{
	for (int i = 0; i <= 200; i++)
	{
		for (int j = 0; j <= 200; j++)
		{
			Matrix4x4 m = Matrix4x4();
			m.rotate(Scene::getRand()*360.f, 0, 1, 0);
			m.scale(Scene::getRand()*0.3+0.85, Scene::getRand()*0.3+0.7, Scene::getRand()*0.3+0.85);
			m.translate(-2+i*(Scene::getRand()*0.2+0.2), 0, 3-j*(Scene::getRand()*0.2+0.2));
			ProxyObject* po = new ProxyObject(o, b, m);
			po->setDisplayNum(n);
			g_scene->addObject(po);
		}
	}
}

void makeTrees(Objects* o, BVH* b, int n)
{
	for (int i = 0; i <= 200; i++)
	{
		float x, z;
		do 
		{
			x = Scene::getRand();
			z = Scene::getRand();
		} while (x*x + z*z > 1.f && sqrtf(x*x + z*z) < 0.125f && x < 0.f && z > 0.f);
		Matrix4x4 m = Matrix4x4();
		m.rotate(Scene::getRand()*360.f, 0, 1, 0);
		m.scale(Scene::getRand()*0.3+0.85, Scene::getRand()*0.3+0.85, Scene::getRand()*0.3+0.85);
		m.translate(x*8.f*100.f, Scene::getRand()*0.5f-0.5f, -z*8.f*100.f);
		if (m.m14 < 100.f && m.m34 > -100.f)
		{
			continue;
		}
		ProxyObject* po = new ProxyObject(o, b, m);
		po->setDisplayNum(n);
		g_scene->addObject(po);
	}
}

void makeFlowers(Objects* o, BVH* b, int n)
{
	for (int i = 0; i <= 390; i++)
	{
		float x, z;
		do 
		{
			x = Scene::getRand();
			z = Scene::getRand();
		} while (x*x + z*z > 1.f);
		Matrix4x4 m = Matrix4x4();
		m.rotate(Scene::getRand()*360.f, 0, 1, 0);
		m *= Matrix4x4().rotateX(Scene::getRand()*20.f+10.f);
		m.scale(Scene::getRand()*0.2+0.9, Scene::getRand()*0.2+0.95, Scene::getRand()*0.2+0.9);
		m.translate(g_camera->eye().x+x*10.f, Scene::getRand()*0.05f-0.025f, g_camera->eye().z-z*10.f);
		ProxyObject* po = new ProxyObject(o, b, m);
		po->setDisplayNum(n);
		g_scene->addObject(po);
	}
}

void addProxyObj(Objects* o, BVH* b, Matrix4x4 &m, int n)
{
	ProxyObject* po = new ProxyObject(o, b, m);
	po->setDisplayNum(n);
	g_scene->addObject(po);
}

void camera01Settings()
{
	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-1.277, 0.158, 2.139));
	g_camera->setLookAt(0.294, 0.511, 0.503);
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(39);
	g_camera->m_aperture = 0.0018f;
	g_camera->m_focusPlane = 2.f;
	g_camera->setShutterSpeed(0.1f);
}

void camera02Settings()
{
	// set up the camera
	g_scene->setBGColor(Vector3(0));
	g_camera->setEye(Vector3(-0.386, 0.372, 1.238));
	g_camera->setLookAt(0.294, 0.511, 0.503);
	g_camera->setUp(Vector3(0, 1, 0));
	g_camera->setFOV(39);
	g_camera->m_aperture = 0.0025f;
	g_camera->m_focusPlane = 0.7f;
	g_camera->setShutterSpeed(0.005f);
}

void makeFinalScene()
{
	g_camera = (Camera*)_aligned_malloc(sizeof(Camera), 16);
	g_scene = new Scene;
	g_image = new Image;

	g_image->resize(1920, 1080);

	g_scene->m_pathTrace = false;
	g_scene->m_numPaths = 1;
	g_scene->m_maxBounces = 5;
	g_scene->m_minSubdivs = 3;
	g_scene->m_maxSubdivs = 5;
	g_scene->setNoise(0.01f);

	camera01Settings();

	//make a raw image from hdr file
	RawImage* BGImage = new RawImage();
	BGImage->loadImage("Textures/hdrvfx_nyany_1_n2_v101_Bg.tga");
	Texture* BGTex = new Texture(BGImage);
	g_scene->setEnvMap(BGTex);
	g_scene->setEnvExposure(1.5f);

	//make a raw image from hdr file
	RawImage* hdrImage2 = new RawImage();
	hdrImage2->loadImage("Images/sky.hdr");
	Texture* hdrTex2 = new Texture(hdrImage2);

	DomeLight* domeLight = new DomeLight; 
	domeLight->setTexture(hdrTex2);
	domeLight->setPower(0.15f);
	domeLight->setSamples(6);
	g_scene->addLight(domeLight);

	Blinn* mat = new Blinn(Vector3(0.9, 0.9, 0.9));
	mat->setSpecExp(30.0f);
	mat->setSpecAmt(0);
	mat->setIor(1.56f);
	mat->setReflectAmt(1.0f);
	mat->setRefractAmt(1.0f);
	mat->setReflectGloss(1.f);
	mat->m_disperse = true;

	RawImage* grassImg = new RawImage();
	grassImg->loadImage("Textures/grassblade2.tga");
	Texture* grassTex = new Texture(grassImg);

	Blinn* grassMtl = new Blinn(Vector3(0.5));
	grassMtl->setSpecExp(20.0f);
	grassMtl->setSpecAmt(0.8);
	grassMtl->setColorMap(grassTex);

	RawImage* dirtImg = new RawImage();
	dirtImg->loadImage("Textures/ground-dirt-texture.tga");
	Texture* dirtTex = new Texture(dirtImg);

	Blinn* dirtMat = new Blinn(Vector3(0.1, 0.1, 0.1));
	dirtMat->setSpecExp(30.0f);
	dirtMat->setSpecAmt(0);
	dirtMat->setIor(1.80f);
	dirtMat->setReflectAmt(0.0f);
	dirtMat->setRefractAmt(0.0f);
	dirtMat->setReflectGloss(1.f);
	dirtMat->setColorMap(dirtTex);

	TriangleMesh *mesh = new TriangleMesh;
	TriangleMesh *mesh2 = new TriangleMesh;
	mesh->load("Models/Final/explosion01.obj");
	mesh2->load("Models/Final/explosion02.obj");
	makeMBMeshObjs(mesh, mesh2, mat);

	RawImage* bulletImg = new RawImage();
	bulletImg->loadImage("Textures/bw2.tga");
	Texture* bulletTex = new Texture(bulletImg);

	Blinn* cBallMat = new Blinn(Vector3(0.01, 0.01, 0.01));
	cBallMat->setSpecExp(15.0f);
	cBallMat->setSpecAmt(0.5);
	cBallMat->setIor(1.80f);
	cBallMat->setReflectAmt(0.0f);
	cBallMat->setRefractAmt(0.0f);
	cBallMat->setReflectGloss(0.9f);
	cBallMat->setColorMap(bulletTex);

	TriangleMesh *cmesh = new TriangleMesh;
	TriangleMesh *cmesh2 = new TriangleMesh;
	cmesh->load("Models/Final/cannonBallT1.obj");
	cmesh2->load("Models/Final/cannonBallT2.obj");
	makeMBMeshObjs(cmesh, cmesh2, cBallMat);

	TriangleMesh *groundPlane = new TriangleMesh;
	groundPlane->load("Models/Final/groundPlane.obj");
	makeMeshObjs(groundPlane, dirtMat);

	//////////////////////////////////////////////////////////////////////////

	Matrix4x4 T02m01 = Matrix4x4();
	T02m01.rotate(0, 0, 1, 0);
	T02m01.scale(0.64,0.64,0.64);
	T02m01.translate(62.872, 0, -27.025);

	RawImage* Tree02BodyImg = new RawImage();
	Tree02BodyImg->loadImage("Textures/AL04brk.tga");
	Texture* Tree02BodyTex = new Texture(Tree02BodyImg);

	Blinn* Tree02BodyMtl = new Blinn(Vector3(0.5));
	Tree02BodyMtl->setSpecExp(20.0f);
	Tree02BodyMtl->setSpecAmt(0.8);
	Tree02BodyMtl->setColorMap(Tree02BodyTex);

	TriangleMesh *Tree02Body = new TriangleMesh;
	Tree02Body->load("Models/Final/Tree02Body.obj");

	RawImage* tree02LeavesImg = new RawImage();
	tree02LeavesImg->loadImage("Textures/AL04aut.tga");
	Texture* tree02LeavesTex = new Texture(tree02LeavesImg);

	Blinn* tree02LeavesMtl = new Blinn(Vector3(0.5));
	tree02LeavesMtl->setSpecExp(20.0f);
	tree02LeavesMtl->setSpecAmt(0.8);
	tree02LeavesMtl->setTranslucency(0.6f);
	tree02LeavesMtl->setColorMap(tree02LeavesTex);
	tree02LeavesMtl->setAlphaMap(tree02LeavesTex);

	TriangleMesh *tree02Leaves = new TriangleMesh;
	tree02Leaves->load("Models/Final/tree02Leaves.obj");
	
	BVH* tree02BVH = new BVH;
	Objects* tree02Os = new Objects;

	Material** tree02Mtl = new Material*[2];
	tree02Mtl[0] = Tree02BodyMtl;
	tree02Mtl[1] = tree02LeavesMtl;

	TriangleMesh** tree02Mesh = new TriangleMesh*[2];
	tree02Mesh[0] = Tree02Body;
	tree02Mesh[1] = tree02Leaves;

	ProxyObject::setupMultiProxy(tree02Mesh, 2, tree02Mtl, tree02Os, tree02BVH);

	makeTrees(tree02Os, tree02BVH, 200);
	addProxyObj(tree02Os, tree02BVH, T02m01, 1500);

	//////////////////////////////////////////////////////////////////////////

	Matrix4x4 T01m01 = Matrix4x4();
	T01m01.rotate(0, 0, 1, 0);
	T01m01.scale(1,1,1);
	T01m01.translate(0, 0, -21.013);

	Matrix4x4 T01m02 = Matrix4x4();
	T01m02.rotate(-105.05, 0, 1, 0);
	T01m02.scale(1, 1, 1);
	T01m02.translate(43.078, 0, -9.234);

	Matrix4x4 T01m03 = Matrix4x4();
	T01m03.rotate(-173.91, 0, 1, 0);
	T01m03.scale(1.164, 1.164, 1.164);
	T01m03.translate(93.86, 0, -53.41);

	Matrix4x4 T01m04 = Matrix4x4();
	T01m04.rotate(100, 0, 1, 0);
	T01m04.scale(0.71, 0.71, 0.71);
	T01m04.translate(10.92, 0, -53.16);

	RawImage* tree01LeavesImg = new RawImage();
	tree01LeavesImg->loadImage("Textures/ML16lef1.tga");
	Texture* tree01LeavesTex = new Texture(tree01LeavesImg);

	Blinn* tree01LeavesMtl = new Blinn(Vector3(0.5));
	tree01LeavesMtl->setSpecExp(20.0f);
	tree01LeavesMtl->setSpecAmt(0.8);
	tree01LeavesMtl->setTranslucency(0.6f);
	tree01LeavesMtl->setColorMap(tree01LeavesTex);
	tree01LeavesMtl->setAlphaMap(tree01LeavesTex);

	TriangleMesh *tree01Leaves = new TriangleMesh;
	tree01Leaves->load("Models/Final/tree01Leaves.obj");

	RawImage* Tree01BodyImg = new RawImage();
	Tree01BodyImg->loadImage("Textures/ML16brk.tga");
	Texture* Tree01BodyTex = new Texture(Tree01BodyImg);

	Blinn* Tree01BodyMtl = new Blinn(Vector3(0.5));
	Tree01BodyMtl->setSpecExp(20.0f);
	Tree01BodyMtl->setSpecAmt(0.8);
	Tree01BodyMtl->setColorMap(Tree01BodyTex);

	TriangleMesh *Tree01Body = new TriangleMesh;
	Tree01Body->load("Models/Final/Tree01Body.obj");

	BVH* tree01BVH = new BVH;
	Objects* tree01Os = new Objects;

	Material** tree01Mtl = new Material*[2];
	tree01Mtl[0] = Tree01BodyMtl;
	tree01Mtl[1] = tree01LeavesMtl;

	TriangleMesh** tree01Mesh = new TriangleMesh*[2];
	tree01Mesh[0] = Tree01Body;
	tree01Mesh[1] = tree01Leaves;

	ProxyObject::setupMultiProxy(tree01Mesh, 2, tree01Mtl, tree01Os, tree01BVH);

	makeTrees(tree01Os, tree01BVH, 1500);
	makeTrees(tree01Os, tree01BVH, 1500);
	addProxyObj(tree01Os, tree01BVH, T01m01, 200);
	addProxyObj(tree01Os, tree01BVH, T01m02, 200);
	addProxyObj(tree01Os, tree01BVH, T01m03, 200);
	addProxyObj(tree01Os, tree01BVH, T01m04, 200);

	//////////////////////////////////////////////////////////////////////////

	TriangleMesh *Tree04Body = new TriangleMesh;
	Tree04Body->load("Models/Final/Tree04Body.obj");

	TriangleMesh *tree04Leaves = new TriangleMesh;
	tree04Leaves->load("Models/Final/tree04Leaves.obj");

	BVH* tree04BVH = new BVH;
	Objects* tree04Os = new Objects;

	TriangleMesh** tree04Mesh = new TriangleMesh*[2];
	tree04Mesh[0] = Tree04Body;
	tree04Mesh[1] = tree04Leaves;

	ProxyObject::setupMultiProxy(tree04Mesh, 2, tree01Mtl, tree04Os, tree04BVH);

	addProxyObj(tree04Os, tree04BVH, Matrix4x4(), 200);

	//////////////////////////////////////////////////////////////////////////

	RawImage* Tree03BodyImg = new RawImage();
	Tree03BodyImg->loadImage("Textures/AL17brk.tga");
	Texture* Tree03BodyTex = new Texture(Tree03BodyImg);

	Blinn* Tree03BodyMtl = new Blinn(Vector3(0.5));
	Tree03BodyMtl->setSpecExp(20.0f);
	Tree03BodyMtl->setSpecAmt(0.8);
	Tree03BodyMtl->setColorMap(Tree03BodyTex);

	TriangleMesh *Tree03Body = new TriangleMesh;
	Tree03Body->load("Models/Final/Tree03Body.obj");
	makeMeshObjs(Tree03Body, Tree03BodyMtl);

	RawImage* tree03LeavesImg = new RawImage();
	tree03LeavesImg->loadImage("Textures/AL17aut.tga");
	Texture* tree03LeavesTex = new Texture(tree03LeavesImg);

	Blinn* tree03LeavesMtl = new Blinn(Vector3(0.5));
	tree03LeavesMtl->setSpecExp(20.0f);
	tree03LeavesMtl->setSpecAmt(0.8);
	tree03LeavesMtl->setTranslucency(0.6f);
	tree03LeavesMtl->setColorMap(tree03LeavesTex);
	tree03LeavesMtl->setAlphaMap(tree03LeavesTex);

	TriangleMesh *tree03Leaves = new TriangleMesh;
	tree03Leaves->load("Models/Final/tree03Leaves.obj");
	makeMeshObjs(tree03Leaves, tree03LeavesMtl);

	//////////////////////////////////////////////////////////////////////////

	Matrix4x4 rotateM = Matrix4x4();
	Matrix4x4 fl02m01 = Matrix4x4();
	fl02m01 = rotateM.rotateZ(5.71) * rotateM.rotateY(90.472) * rotateM.rotateX(27.652);
	fl02m01.translate(-1.139, 0.013, 1.801);

	RawImage* flower02BodyImg = new RawImage();
	flower02BodyImg->loadImage("Textures/grass-color-23.tga");
	Texture* flower02BodyTex = new Texture(flower02BodyImg);

	Blinn* flower02BodyMtl = new Blinn(Vector3(0.5));
	flower02BodyMtl->setSpecExp(10.0f);
	flower02BodyMtl->setSpecAmt(0.5);
	flower02BodyMtl->setColorMap(flower02BodyTex);

	TriangleMesh *flower02Body = new TriangleMesh;
	flower02Body->load("Models/Final/flower02Body.obj");

	RawImage* flower02BulbImg = new RawImage();
	flower02BulbImg->loadImage("Textures/bud-yellow-1.tga");
	Texture* flower02BulbTex = new Texture(flower02BulbImg);

	RawImage* flower02BulbNormalImg = new RawImage();
	flower02BulbNormalImg->loadImage("Textures/bud-yellow-1-bump_NRM.tga");
	Texture* flower02BulbNormalTex = new Texture(flower02BulbNormalImg);

	Blinn* flower02BulbMtl = new Blinn(Vector3(0.5));
	flower02BulbMtl->setSpecExp(1.0f);
	flower02BulbMtl->setSpecAmt(0.0);
	flower02BulbMtl->setColorMap(flower02BulbTex);
	flower02BulbMtl->setNormalMap(flower02BulbNormalTex);

	TriangleMesh *flower02Bulb = new TriangleMesh;
	flower02Bulb->load("Models/Final/flower02Bulb.obj");

	RawImage* flower02LeavesImg = new RawImage();
	flower02LeavesImg->loadImage("Textures/grass-color-18.tga");
	Texture* flower02LeavesTex = new Texture(flower02LeavesImg);

	Blinn* flower02LeavesMtl = new Blinn(Vector3(0.5));
	flower02LeavesMtl->setSpecExp(20.0f);
	flower02LeavesMtl->setSpecAmt(0.5);
	flower02LeavesMtl->setTranslucency(0.5);
	flower02LeavesMtl->setColorMap(flower02LeavesTex);

	TriangleMesh *flower02Leaves = new TriangleMesh;
	flower02Leaves->load("Models/Final/flower02Leaves.obj");

	RawImage* flower02PinkPetalsImg = new RawImage();
	flower02PinkPetalsImg->loadImage("Textures/petal-pink-02.tga");
	Texture* flower02PinkPetalsTex = new Texture(flower02PinkPetalsImg);

	Blinn* flower02PinkPetalsMtl = new Blinn(Vector3(0.5));
	flower02PinkPetalsMtl->setSpecExp(10.0f);
	flower02PinkPetalsMtl->setSpecAmt(0.3);
	flower02PinkPetalsMtl->setTranslucency(0.6);
	flower02PinkPetalsMtl->setColorMap(flower02PinkPetalsTex);

	RawImage* flower02YellowPetalsImg = new RawImage();
	flower02YellowPetalsImg->loadImage("Textures/petal-yellow-1.tga");
	Texture* flower02YellowPetalsTex = new Texture(flower02YellowPetalsImg);

	Blinn* flower02YellowPetalsMtl = new Blinn(Vector3(0.5));
	flower02YellowPetalsMtl->setSpecExp(10.0f);
	flower02YellowPetalsMtl->setSpecAmt(0.3);
	flower02YellowPetalsMtl->setTranslucency(0.6);
	flower02YellowPetalsMtl->setColorMap(flower02YellowPetalsTex);

	RawImage* flower02WhitePetalsImg = new RawImage();
	flower02WhitePetalsImg->loadImage("Textures/petal-white-3.tga");
	Texture* flower02WhitePetalsTex = new Texture(flower02WhitePetalsImg);

	Blinn* flower02WhitePetalsMtl = new Blinn(Vector3(0.5));
	flower02WhitePetalsMtl->setSpecExp(10.0f);
	flower02WhitePetalsMtl->setSpecAmt(0.3);
	flower02WhitePetalsMtl->setTranslucency(0.6);
	flower02WhitePetalsMtl->setColorMap(flower02WhitePetalsTex);

	TriangleMesh *flower02Petals = new TriangleMesh;
	flower02Petals->load("Models/Final/flower02Petals.obj");

	BVH* flower02PinkBVH = new BVH;
	Objects* flower02PinkOs = new Objects;
	BVH* flower02YellowBVH = new BVH;
	Objects* flower02YellowOs = new Objects;
	BVH* flower02WhiteBVH = new BVH;
	Objects* flower02WhiteOs = new Objects;

	Material** flower02PinkMtl = new Material*[4];
	flower02PinkMtl[3] = flower02BodyMtl;
	flower02PinkMtl[2] = flower02BulbMtl;
	flower02PinkMtl[1] = flower02LeavesMtl;
	flower02PinkMtl[0] = flower02PinkPetalsMtl;

	Material** flower02YellowMtl = new Material*[4];
	flower02YellowMtl[3] = flower02BodyMtl;
	flower02YellowMtl[2] = flower02BulbMtl;
	flower02YellowMtl[1] = flower02LeavesMtl;
	flower02YellowMtl[0] = flower02YellowPetalsMtl;

	Material** flower02WhiteMtl = new Material*[4];
	flower02WhiteMtl[3] = flower02BodyMtl;
	flower02WhiteMtl[2] = flower02BulbMtl;
	flower02WhiteMtl[1] = flower02LeavesMtl;
	flower02WhiteMtl[0] = flower02WhitePetalsMtl;

	TriangleMesh** flower02Mesh = new TriangleMesh*[4];
	flower02Mesh[3] = flower02Body;
	flower02Mesh[2] = flower02Bulb;
	flower02Mesh[1] = flower02Leaves;
	flower02Mesh[0] = flower02Petals;

	ProxyObject::setupMultiProxy(flower02Mesh, 4, flower02PinkMtl, flower02PinkOs, flower02PinkBVH);
	ProxyObject::setupMultiProxy(flower02Mesh, 4, flower02YellowMtl, flower02YellowOs, flower02YellowBVH);
	ProxyObject::setupMultiProxy(flower02Mesh, 4, flower02WhiteMtl, flower02WhiteOs, flower02WhiteBVH);

	addProxyObj(flower02PinkOs, flower02PinkBVH, fl02m01, 10);
	makeFlowers(flower02PinkOs, flower02PinkBVH, 2);
	makeFlowers(flower02YellowOs, flower02YellowBVH, 2);
	makeFlowers(flower02WhiteOs, flower02WhiteBVH, 2);

	//////////////////////////////////////////////////////////////////////////

	Matrix4x4 fl01m01 = Matrix4x4();
	fl01m01.rotate(0, 0, 1, 0);
	fl01m01.scale(0.6703, 0.6703, 0.6703);
	fl01m01.translate(-1.014, 0, 1.302);

	Matrix4x4 fl01m02 = Matrix4x4();
	fl01m02.rotate(-87.07, 0, 1, 0);
	fl01m02.scale(0.54, 0.54, 0.54);
	fl01m02.translate(-0.464, 0, 0.149);

	Matrix4x4 fl01m03 = Matrix4x4();
	fl01m03.rotate(0, 0, 1, 0);
	fl01m03.scale(0.88487, 0.88487, 0.88487);
	fl01m03.translate(1.264, 0, 0.207);

	Matrix4x4 fl01m04 = Matrix4x4();
	fl01m04.rotate(0, 0, 1, 0);
	fl01m04.scale(0.67, 0.67, 0.67);
	fl01m04.translate(1.96, 0, 1.009);

	RawImage* flower01BigLeavesImg = new RawImage();
	flower01BigLeavesImg->loadImage("Textures/FL30lef1.tga");
	Texture* flower01BigLeavesTex = new Texture(flower01BigLeavesImg);

	Blinn* flower01BigLeavesMtl = new Blinn(Vector3(0.5));
	flower01BigLeavesMtl->setSpecExp(20.0f);
	flower01BigLeavesMtl->setSpecAmt(0.8);
	flower01BigLeavesMtl->setTranslucency(0.6f);
	flower01BigLeavesMtl->setColorMap(flower01BigLeavesTex);
	flower01BigLeavesMtl->setAlphaMap(flower01BigLeavesTex);

	TriangleMesh *flower01BigLeaves = new TriangleMesh;
	flower01BigLeaves->load("Models/Final/flower01BigLeaves.obj");

	RawImage* flower01BodyImg = new RawImage();
	flower01BodyImg->loadImage("Textures/FL30stm1.tga");
	Texture* flower01BodyTex = new Texture(flower01BodyImg);

	Blinn* flower01BodyMtl = new Blinn(Vector3(0.5));
	flower01BodyMtl->setSpecExp(20.0f);
	flower01BodyMtl->setSpecAmt(0.8);
	flower01BodyMtl->setColorMap(flower01BodyTex);

	TriangleMesh *flower01Body = new TriangleMesh;
	flower01Body->load("Models/Final/flower01Body.obj");

	RawImage* flower01Bulbs01Img = new RawImage();
	flower01Bulbs01Img->loadImage("Textures/FL30flo1.tga");
	Texture* flower01Bulbs01Tex = new Texture(flower01Bulbs01Img);

	Blinn* flower01Bulbs01Mtl = new Blinn(Vector3(0.5));
	flower01Bulbs01Mtl->setSpecExp(20.0f);
	flower01Bulbs01Mtl->setSpecAmt(0.8);
	flower01Bulbs01Mtl->setColorMap(flower01Bulbs01Tex);

	TriangleMesh *flower01Bulbs01 = new TriangleMesh;
	flower01Bulbs01->load("Models/Final/flower01Bulbs01.obj");

	Blinn* flower01Bulbs02Mtl = new Blinn(Vector3(0.5));
	flower01Bulbs02Mtl->setSpecExp(20.0f);
	flower01Bulbs02Mtl->setSpecAmt(0.8);
	flower01Bulbs02Mtl->setColorMap(flower01BodyTex);

	TriangleMesh *flower01Bulbs02 = new TriangleMesh;
	flower01Bulbs02->load("Models/Final/flower01Bulbs02.obj");

	Blinn* flower01Bulbs03Mtl = new Blinn(Vector3(1.0, 0.64, 0.15));
	flower01Bulbs03Mtl->setSpecExp(20.0f);
	flower01Bulbs03Mtl->setSpecAmt(0.8);

	TriangleMesh *flower01Bulbs03 = new TriangleMesh;
	flower01Bulbs03->load("Models/Final/flower01Bulbs03.obj");

	RawImage* flower01PetalsImg = new RawImage();
	flower01PetalsImg->loadImage("Textures/FL30pet1.tga");
	Texture* flower01PetalsTex = new Texture(flower01PetalsImg);

	Blinn* flower01PetalsMtl = new Blinn(Vector3(0.5));
	flower01PetalsMtl->setSpecExp(20.0f);
	flower01PetalsMtl->setSpecAmt(0.8);
	flower01PetalsMtl->setTranslucency(0.6f);
	flower01PetalsMtl->setColorMap(flower01PetalsTex);

	TriangleMesh *flower01Petals = new TriangleMesh;
	flower01Petals->load("Models/Final/flower01Petals.obj");

	RawImage* flower01PistilsImg = new RawImage();
	flower01PistilsImg->loadImage("Textures/FL30stm2.tga");
	Texture* flower01PistilsTex = new Texture(flower01PistilsImg);

	Blinn* flower01PistilsMtl = new Blinn(Vector3(0.5));
	flower01PistilsMtl->setSpecExp(20.0f);
	flower01PistilsMtl->setSpecAmt(0.8);
	flower01PistilsMtl->setColorMap(flower01PistilsTex);

	TriangleMesh *flower01Pistils = new TriangleMesh;
	flower01Pistils->load("Models/Final/flower01Pistils.obj");

	RawImage* flower01SmallLeavesImg = new RawImage();
	flower01SmallLeavesImg->loadImage("Textures/FL30lef2.tga");
	Texture* flower01SmallLeavesTex = new Texture(flower01SmallLeavesImg);

	Blinn* flower01SmallLeavesMtl = new Blinn(Vector3(0.5));
	flower01SmallLeavesMtl->setSpecExp(20.0f);
	flower01SmallLeavesMtl->setSpecAmt(0.8);
	flower01SmallLeavesMtl->setTranslucency(0.6f);
	flower01SmallLeavesMtl->setColorMap(flower01SmallLeavesTex);
	flower01SmallLeavesMtl->setAlphaMap(flower01SmallLeavesTex);

	TriangleMesh *flower01SmallLeaves = new TriangleMesh;
	flower01SmallLeaves->load("Models/Final/flower01SmallLeaves.obj");

	BVH* flower01BVH = new BVH;
	Objects* flower01Os = new Objects;

	Material** flower01Mtl = new Material*[8];
	flower01Mtl[0] = flower01BigLeavesMtl;
	flower01Mtl[1] = flower01BodyMtl;
	flower01Mtl[2] = flower01Bulbs01Mtl;
	flower01Mtl[3] = flower01Bulbs02Mtl;
	flower01Mtl[4] = flower01Bulbs03Mtl;
	flower01Mtl[5] = flower01PetalsMtl;
	flower01Mtl[6] = flower01PistilsMtl;
	flower01Mtl[7] = flower01SmallLeavesMtl;

	TriangleMesh** flower01Mesh = new TriangleMesh*[8];
	flower01Mesh[0] = flower01BigLeaves;
	flower01Mesh[1] = flower01Body;
	flower01Mesh[2] = flower01Bulbs01;
	flower01Mesh[3] = flower01Bulbs02;
	flower01Mesh[4] = flower01Bulbs03;
	flower01Mesh[5] = flower01Petals;
	flower01Mesh[6] = flower01Pistils;
	flower01Mesh[7] = flower01SmallLeaves;

	ProxyObject::setupMultiProxy(flower01Mesh, 8, flower01Mtl, flower01Os, flower01BVH);

	addProxyObj(flower01Os, flower01BVH, fl01m01, 10);
	addProxyObj(flower01Os, flower01BVH, fl01m02, 10);
	addProxyObj(flower01Os, flower01BVH, fl01m03, 10);
	addProxyObj(flower01Os, flower01BVH, fl01m04, 10);

	//////////////////////////////////////////////////////////////////////////

	TriangleMesh *grass = new TriangleMesh;
	grass->load("Models/testGrass2.obj");
	BVH* b = new BVH;
	Objects* o = new Objects;
	ProxyObject::setupProxy(grass, grassMtl, o, b);

	makeProxyGrid(o, b, 5000);

	// let objects do pre-calculations if needed
	g_scene->preCalc();
}