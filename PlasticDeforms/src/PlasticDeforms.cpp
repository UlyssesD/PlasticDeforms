#include <stdlib.h>
#include <iostream>
#include <cstdio>
using namespace std;

//VegaFEM inclusions
#include "openGL-headers.h"
#include "lighting.h"
#include "objMesh.h"
#include "sceneObjectDeformable.h"

//Project inclusions
#include "initGraphics.h"


//Program variables
string meshFilename = "";
ObjMesh* mesh = NULL;
Lighting* lighting = NULL;
SceneObjectDeformable* sceneObj = NULL;

string lightingConfigFilename = "..\\models\\cube\\cube.lighting";
char windowTitleBase[4096] = "Plastic Deformations - Toy Example";
void displayFunction(void);
int windowID;
int windowWidth = 800;
int windowHeight = 600;
double zNear = 0.01;
double zFar = 10.0;
double cameraRadius = 5;
double focusPositionX, focusPositionY, focusPositionZ = 0.0;
double cameraLongitude = 45.0, cameraLatitude = 45.0;
SphericalCamera * camera = NULL;


int g_iMiddleMouseButton = 0, g_iRightMouseButton = 0;
int g_vMousePos[2] = { 0,0 };
bool renderWireframe = false, renderNormals = false;

// graphics loop function.
void displayFunction(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Setup model transformations.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	camera->Look();
	
	sceneObj->SetLighting(lighting);

	lighting->LightScene();

	//render the current mesh
	glEnable(GL_LIGHTING);
	sceneObj->Render();

	glDisable(GL_LIGHTING);

	//if enabled, render wireframe
	if(renderWireframe)
		sceneObj->RenderEdges();

	//if enabled, render face normals
	if (renderNormals)
		sceneObj->RenderNormals();

	glutSwapBuffers();
}

void mouseButtonActivityFunction(int button, int state, int x, int y)
{
	switch (button)
	{
		case GLUT_MIDDLE_BUTTON:
			g_iMiddleMouseButton = (state == GLUT_DOWN);
			break;

		case GLUT_RIGHT_BUTTON:
			g_iRightMouseButton = (state == GLUT_DOWN);
			break;
	}

	g_vMousePos[0] = x;
	g_vMousePos[1] = y;
}

void mouseMotionFunction(int x, int y)
{
	int mouseDeltaX = x - g_vMousePos[0];
	int mouseDeltaY = y - g_vMousePos[1];

	g_vMousePos[0] = x;
	g_vMousePos[1] = y;

	if (g_iRightMouseButton) // right mouse button handles camera rotations
	{
		const double factor = 0.2;
		camera->MoveRight(factor * mouseDeltaX);
		camera->MoveUp(factor * mouseDeltaY);
	}

	if (g_iMiddleMouseButton) // handle zoom in/out
	{
		const double factor = 0.1;
		camera->ZoomIn(cameraRadius * factor * mouseDeltaY);
	}
}

void idleFunction() 
{
	glutPostRedisplay();
}

void reshape(int x, int y)
{
	glViewport(0, 0, x, y);

	windowWidth = x;
	windowHeight = y;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// compute the window aspect ratio 
	gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, zNear, zFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void keyboardFunction(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 27:
		case 'q':
			exit(0);
			break;

		case 'w':
			renderWireframe = !renderWireframe;
			break;

		case 'n':
			renderNormals = !renderNormals;
			break;

		case '\\':
			camera->Reset();
			break;
	}
}

int main(int argc, char* argv[]) 
{

	if (argc < 2)
	{
		printf("Missing input mesh .obj \n");
		printf("Usage: %s [.obj file]. \n", argv[0]);
		return 1;
	}

	printf("Starting application. \n");

	meshFilename = string(argv[1]);
	printf("Loading file %s \n", meshFilename.c_str());

	//load mesh .obj file
	mesh = new ObjMesh(meshFilename);

	if (mesh == NULL)
		printf("Error: failed to load input mesh. \n");
	else
		printf("Success: Number of vertices: %d, Number of faces: %d \n", mesh->getNumVertices(), mesh->getNumFaces());
	
	mesh->buildFaceNormals();

	sceneObj = new SceneObjectDeformable(mesh);

	initGLUT(argc, argv, windowTitleBase, windowWidth, windowHeight, &windowID);
	initGraphics(windowWidth, windowHeight); // more OpenGL initialization calls
	
	// init lighting
	try
	{
		lighting = new Lighting(lightingConfigFilename.c_str());
	}
	catch (int exceptionCode)
	{
		printf("Error (%d) reading lighting information from %s .\n", exceptionCode, lightingConfigFilename.c_str());
		exit(1);
	}

	// init camera
	delete(camera);
	double virtualToPhysicalPositionFactor = 1.0;
	initCamera(cameraRadius, cameraLongitude, cameraLatitude, focusPositionX, focusPositionY, focusPositionZ, 1.0 / virtualToPhysicalPositionFactor, &zNear, &zFar, &camera);

	//try to visualize the mesh
	//meshRender->render(OBJMESHRENDER_TRIANGLES, OBJMESHRENDER_MATERIAL);
	glutMainLoop(); // you have reached the point of no return..

	return 0;
}