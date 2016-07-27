#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <limits>
#include <direct.h>
using namespace std;

//VegaFEM inclusions
#include "openGL-headers.h"
#include "openGLHelper.h"
#include "GL/glui.h"
#include "configFile.h"
#include "lighting.h"
#include "matrixIO.h"
#include "modalMatrix.h"
#include "sceneObjectReducedCPU.h"
#include "reducedStVKForceModel.h"
#include "implicitNewmarkDense.h"

//Project inclusions
#include "initGraphics.h"


//funcion declaration
void displayFunction(void);

//Program variables
string configFilesDir;
string configFilename, outputFilename;
string extraSceneGeometryFilename;

ObjMesh * mesh = NULL;
Lighting * lighting = NULL;
SceneObjectDeformable * extraSceneGeometry = NULL;
SceneObjectReduced * deformableObjectRenderingMeshReduced = NULL;
SceneObjectReducedCPU* deformableObjectRenderingMeshCPU = NULL;
ModalMatrix* renderingModalMatrix = NULL;
ImplicitNewmarkDense * implicitNewmarkDense = NULL;
StVKReducedInternalForces * stVKReducedInternalForces = NULL;
StVKReducedStiffnessMatrix * stVKReducedStiffnessMatrix = NULL;
ReducedForceModel * reducedForceModel = NULL;
ReducedStVKForceModel * reducedStVKForceModel;
GLUI * glui;

char windowTitleBase[4096] = "Plastic Deformations - Toy Example";
void displayFunction(void);
int windowID;
int windowWidth = 800;
int windowHeight = 600;
double zNear = 0.01;
double zFar = 10.0;
double cameraRadius = 5;
double focusPositionX = 0, focusPositionY = 0, focusPositionZ = 0;
double cameraLongitude, cameraLattitude;
//double threshold = numeric_limits<double>::infinity();
double threshold = 1E-4;
SphericalCamera * camera = NULL;


int g_iMiddleMouseButton = 0, g_iRightMouseButton = 0;
int g_vMousePos[2] = { 0,0 };


// options for the config file
char deformableObjectFilename[4096];
char modesFilename[4096];
char cubicPolynomialFilename[4096];
char lightingConfigFilename[4096];
float dampingMassCoef = 0.0;
float dampingStiffnessCoef = 0.0;
char backgroundColorString[4096] = "255 255 255";
int renderOnGPU;
int substepsPerTimeStep = 1;
int plasticDeformationsEnabled = 1;
float plasticThreshold = 1E9;
float frequencyScaling = 1;

int renderWireframe = 0, renderNormals = 0;

int nRendering;
int r;
double * q = NULL;
double * fq = NULL;
double * fqBase = NULL;

float timeStep = 1.0 / 30;
float newmarkBeta = 0.25;
float newmarkGamma = 0.5;

double impulse = 1*pow(10, 6), max_impulse = 100.0, step = 10;

int numverts = 0;
int* vertices;
double* displacements;
int vertex = 5956;

Vec3d trajectory[4] = { Vec3d(-0.03, 0.03, 0), Vec3d(-0.03, 0.015, 0), Vec3d(0.03, 0.015, 0), Vec3d(0.03, 0.03, 0) };

void applyImpulseForce()
{
	printf("Applying an impulse force of %f Newton along the Y axes on selected vertices. \n", impulse);
	
	double externalForce[3] = { 0, -1*impulse, 0 };

	//printf("Externalforce: [%f, %f, %f]\n", externalForce[0], externalForce[1], externalForce[2]);

	double cos = acos(dot(Vec3d(0, impulse, 0), Vec3d(0, 1, 0))/ len(Vec3d(0, impulse, 0)) * len(Vec3d(0, 1, 0)));
	//printf("Is force perpendicular: %d. \n", cos == PI);

	for (int i = 0; i < numverts; i++) {
		//printf("Applying an impulse force of %f Newton along the Y axes on vertex %d. \n", impulse, vertices[i]);
			renderingModalMatrix->ProjectSingleVertex(vertices[i], externalForce[0], externalForce[1], externalForce[2], fq);
			 
			for (int j = 0; j < r; j++)
				fqBase[j] = fqBase[j] + fq[j];
	}
	//renderingModalMatrix->ProjectSingleVertex(vertex,
	//	externalForce[0], externalForce[1], externalForce[2], fq);
	
	memcpy(fq, fqBase, sizeof(double) * r);
	// set the reduced external forces
	implicitNewmarkDense->SetExternalForces(fq);



	//printf("System kinetic energy: %f \n", implicitNewmarkDense->GetKineticEnergy());

	do	{
		
		for (int i = 0; i < substepsPerTimeStep; i++)
		{
			int code = implicitNewmarkDense->DoTimestep();
			if (code != 0)
			{
				printf("The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep.\n");
				implicitNewmarkDense->ResetToRest();
			}
			
			//memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

			// compute u=Uq
			deformableObjectRenderingMeshReduced->Setq(implicitNewmarkDense->Getq());
			deformableObjectRenderingMeshReduced->Compute_uUq();

			deformableObjectRenderingMeshReduced->BuildNormals();

			displayFunction();

		}

	} while (implicitNewmarkDense->GetKineticEnergy() > threshold);

	printf("System kinetic energy: %f \n", implicitNewmarkDense->GetKineticEnergy());
	
	//memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);
	
	// compute u=Uq
	deformableObjectRenderingMeshReduced->Setq(implicitNewmarkDense->Getq());
	deformableObjectRenderingMeshReduced->Compute_uUq();

	//deformableObjectRenderingMeshReduced->BuildNeighboringStructure();
	deformableObjectRenderingMeshReduced->BuildNormals();

	glutPostRedisplay();
	
	for (int i = 0; i < r; i++)
	{
		fq[i] = 0;
		fqBase[i] = 0;
	}
		
}

void exit_buttonCallBack(int code)
{
	exit(0);
}

void impulse_spinnerCallBack(int code)
{
	glui->sync_live();
}

// calls all GLUI callbacks, except the listBox callbacks
void callAllUICallBacks()
{
	impulse_spinnerCallBack(0);
}

void Sync_GLUI()
{
	glui->sync_live();
}

void initGLUI()
{
	glui = GLUI_Master.create_glui("Controls", 0, windowWidth + 10, 0);
	glui->add_spinner("Impulse force:", GLUI_SPINNER_DOUBLE , &impulse, 0, impulse_spinnerCallBack);

	glui->add_button("Exit program", 0, exit_buttonCallBack);

	Sync_GLUI();
	glui->set_main_gfx_window(windowID);
}

// graphics loop function.
void displayFunction(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Setup model transformations.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	camera->Look();
	
	deformableObjectRenderingMeshReduced->SetLighting(lighting);

	lighting->LightScene();

	//render the current mesh
	glEnable(GL_LIGHTING);
	deformableObjectRenderingMeshReduced->Render();
	// render any extra scene geometry
	if (extraSceneGeometry != NULL)
		extraSceneGeometry->Render();

	glDisable(GL_LIGHTING);

	//if enabled, render wireframe
	if(renderWireframe)
		deformableObjectRenderingMeshReduced->RenderEdges();
	
	//if enabled, render face normals
	if (renderNormals)
		deformableObjectRenderingMeshReduced->RenderNormals();

	//DrawArrow(0, 1, 0, 0, 1, 0, 0.05, 0.07);

	glutSwapBuffers();
}

void simulateTool() 
{

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
	Vec3d u; Vec3d  dist;

	switch (key)
	{
		case 27:

		case 'v':
			printf("Calculating impact area;\n");
			numverts = 0;
			vertices = (int *) malloc(deformableObjectRenderingMeshReduced->GetNumVertices() * sizeof(int));
			double v[3];
			for (int i = 0; i < deformableObjectRenderingMeshReduced->GetNumVertices(); i++) {
				deformableObjectRenderingMeshReduced->GetSingleVertexPosition(i, &v[0], &v[1], &v[2]);
				if (v[1] >= 0.0005 && v[0] >= -0.03 && v[0] <= 0.03 && v[2] >= -0.03 && v[2] <= 0.03) {
					//printf("Adding v[%d]: %f, %f, %f; \n", i, v[0], v[1], v[2]);
					vertices[numverts] = i;
					numverts++;
					//renderingModalMatrix->ProjectSingleVertex(i, externalForce[0], externalForce[1], externalForce[2], fq);
				}
			}
			printf("final size of impact area: %d; \n", numverts);
			break;
		case 'q':
			exit(0);
			break;

		case ' ':
			applyImpulseForce();
			break;

		case 'w':
			renderWireframe = !renderWireframe;
			break;

		case 'n':
			renderNormals = !renderNormals;
			break;

		case 'r':
			deformableObjectRenderingMeshReduced->ResetDeformationToRest();
			implicitNewmarkDense->ResetToRest();

			//rebuild normals
			deformableObjectRenderingMeshReduced->BuildNormals();
			if(extraSceneGeometry != NULL) extraSceneGeometry->ResetDeformationToRest();
			break;

		case 's':
			char name[100];
			sprintf(name, "output\\%s_%dN.obj", outputFilename.c_str(), (int) impulse);
			deformableObjectRenderingMeshReduced->GetMesh()->save(name, 1);
			break;
		
		case 't':
			printf("Starting tool movement simulation \n");
			
			printf("moving vertices. \n");

			for (int p = 0; p < size(trajectory) - 1; p++) {
				dist = (trajectory[p + 1] - trajectory[p]); /// len(trajectory[1] - trajectory[0]);

				printf("d: [%f %f %f] \n", dist[0], dist[1], dist[2]);

				displacements = (double *)malloc(3 * extraSceneGeometry->GetNumVertices() * sizeof(double));
				for (int substeps = 0; substeps < 100; substeps++)
				{
					for (int i = 0; i < extraSceneGeometry->GetNumVertices(); i++)
					{
						displacements[i * 3 + 0] = dist[0] / 100;
						displacements[i * 3 + 1] = dist[1] / 100;
						displacements[i * 3 + 2] = dist[2] / 100;
					}
					extraSceneGeometry->AddVertexDeformations(displacements);
					displayFunction();
					simulateTool();
				}

				free(displacements);
			}
			break;

		case 'l':
			_mkdir(("output\\" + outputFilename).c_str());

			//apply impulse force with predefined range
			for (impulse = step; impulse <= max_impulse; impulse = impulse + step)
			{
				//Apply the impulse
				applyImpulseForce();

				//Save the result
				char name[100];
				sprintf(name, "output\\%s\\%s_%dN.obj",outputFilename.c_str(), outputFilename.c_str(), (int) impulse);
				deformableObjectRenderingMeshReduced->GetMesh()->save(name, 1);

				//Reset mesh to rest
				deformableObjectRenderingMeshReduced->ResetDeformationToRest();
				implicitNewmarkDense->ResetToRest();

				//rebuild normals
				//deformableObjectRenderingMeshReduced->BuildNeighboringStructure();
				deformableObjectRenderingMeshReduced->BuildNormals();
			}

			printf("Operation Completed. \n");
			break;

		case '\\':
			camera->Reset();
			break;
	}
}

void initConfigurations()
{
	ConfigFile configFile;

	// specify the entries of the config file
	configFile.addOptionOptional("windowWidth", &windowWidth, 800);
	configFile.addOptionOptional("windowHeight", &windowHeight, 800);

	configFile.addOption("deformableObjectFilename", deformableObjectFilename);
	configFile.addOptionOptional("modesFilename", modesFilename, "__none");
	configFile.addOptionOptional("cubicPolynomialFilename", cubicPolynomialFilename, "__none");

	configFile.addOption("dampingMassCoef", &dampingMassCoef);
	configFile.addOption("dampingStiffnessCoef", &dampingStiffnessCoef);

	configFile.addOptionOptional("plasticThreshold", &plasticThreshold, plasticThreshold);

	//configFile.addOption("deformableObjectCompliance", &deformableObjectCompliance);
	configFile.addOption("frequencyScaling", &frequencyScaling);

	configFile.addOptionOptional("cameraRadius", &cameraRadius, 17.5);
	configFile.addOptionOptional("focusPositionX", &focusPositionX, 0.0);
	configFile.addOptionOptional("focusPositionY", &focusPositionY, 0.0);
	configFile.addOptionOptional("focusPositionZ", &focusPositionZ, 0.0);
	configFile.addOptionOptional("cameraLongitude", &cameraLongitude, 45.0);
	configFile.addOptionOptional("cameraLattitude", &cameraLattitude, 45.0);

	configFile.addOptionOptional("renderWireframe", &renderWireframe, 1);

	//configFile.addOptionOptional("extraSceneGeometry", extraSceneGeometryFilename, "__none");

	configFile.addOptionOptional("backgroundColor", backgroundColorString, backgroundColorString);

	string lightingConfigFilenameDefault = configFilesDir + "default.lighting";
	configFile.addOptionOptional("lightingConfigFilename", lightingConfigFilename,
		(char*)lightingConfigFilenameDefault.c_str());

	configFile.addOptionOptional("substepsPerTimeStep",
		&substepsPerTimeStep, substepsPerTimeStep);

	configFile.addOptionOptional("renderOnGPU", &renderOnGPU, 0);

	// parse the configuration file
	if (configFile.parseOptions((char*)configFilename.c_str()) != 0)
	{
		printf("Error: unable to load the configuration file.\n");
		exit(1);
	}
	// the config variables have now been loaded with their specified values

	// informatively print the variables (with assigned values) that were just parsed
	configFile.printOptions();
}

void initScene()
{
	// init lighting
	try
	{
		lighting = new Lighting(lightingConfigFilename);
	}
	catch (int exceptionCode)
	{
		printf("Error (%d) reading lighting information from %s .\n", exceptionCode, lightingConfigFilename);
		exit(1);
	}

	// init camera
	delete(camera);
	double virtualToPhysicalPositionFactor = 1.0;
	initCamera(cameraRadius, cameraLongitude, cameraLattitude, focusPositionX, focusPositionY, focusPositionZ, 1.0 / virtualToPhysicalPositionFactor, &zNear, &zFar, &camera);

	// load the rendering modes of the deformable object
	int nRendering;
	float * URenderingFloat;
	ReadMatrixFromDisk_(modesFilename, &nRendering, &r, &URenderingFloat);
	nRendering /= 3;
	double * URendering = (double*)malloc(sizeof(double) * 3 * nRendering * r);
	for (int i = 0; i < 3 * nRendering * r; i++)
		URendering[i] = URenderingFloat[i];
	free(URenderingFloat);
	renderingModalMatrix = new ModalMatrix(nRendering, r, URendering);
	free(URendering); // ModalMatrix made an internal copy

	// init room for reduced coordinates and reduced forces
	q = (double*)calloc(r, sizeof(double));
	fq = (double*)calloc(r, sizeof(double));
	fqBase = (double*)calloc(r, sizeof(double));

	// initialize the CPU rendering class for the deformable object
	deformableObjectRenderingMeshCPU = new SceneObjectReducedCPU(deformableObjectFilename, renderingModalMatrix); // uses CPU to compute u=Uq

	deformableObjectRenderingMeshReduced = deformableObjectRenderingMeshCPU;

	deformableObjectRenderingMeshReduced->ResetDeformationToRest();

	// make reduced mass matrix (=identity)
	double * massMatrix = (double*)calloc(r*r, sizeof(double));
	for (int i = 0; i<r; i++)
		massMatrix[ELT(r, i, i)] = 1.0;

	
	stVKReducedInternalForces = new StVKReducedInternalForces(cubicPolynomialFilename); // normal version

	// create stiffness matrix polynomials
	stVKReducedStiffnessMatrix = new StVKReducedStiffnessMatrix(stVKReducedInternalForces);

	// create the "internal force models" 
	reducedStVKForceModel = new ReducedStVKForceModel(stVKReducedInternalForces, stVKReducedStiffnessMatrix);
	//reducedLinearStVKForceModel = new ReducedLinearStVKForceModel(stVKReducedStiffnessMatrix);
	reducedForceModel = reducedStVKForceModel;

	// init the implicit Newmark
	implicitNewmarkDense = new ImplicitNewmarkDense(r, timeStep, massMatrix, reducedForceModel, ImplicitNewmarkDense::positiveDefiniteMatrixSolver, dampingMassCoef, dampingStiffnessCoef);
	implicitNewmarkDense->SetTimestep(timeStep / substepsPerTimeStep);
	implicitNewmarkDense->SetNewmarkBeta(newmarkBeta);
	implicitNewmarkDense->SetNewmarkGamma(newmarkGamma);
	//implicitNewmarkDense->SetInternalForceScalingFactor(frequencyScaling * frequencyScaling);
	deformableObjectRenderingMeshReduced->BuildNeighboringStructure();
	//implicitNewmarkDense->UseStaticSolver(1);
	free(massMatrix);

	if (strcmp(extraSceneGeometryFilename.c_str(),"__none") != 0)
	{
		extraSceneGeometry = new SceneObjectDeformable(extraSceneGeometryFilename.c_str());
		extraSceneGeometry->BuildNeighboringStructure();
		extraSceneGeometry->BuildNormals(85.0);
	}
	else
		extraSceneGeometry = NULL;

	// compute lowest frequency of the system (smallest eigenvalue of K)
	double * K = (double*)malloc(sizeof(double) * r * r);
	double * zero = (double*)calloc(r, sizeof(double));
	stVKReducedStiffnessMatrix->Evaluate(zero, K);

	// set background color
	int colorR, colorG, colorB;
	sscanf(backgroundColorString, "%d %d %d", &colorR, &colorG, &colorB);
	glClearColor(1.0 * colorR / 255, 1.0 * colorG / 255, 1.0 * colorB / 255, 0.0);

	Sync_GLUI();
	callAllUICallBacks();
}

int main(int argc, char* argv[]) 
{

	if (argc < 2)
	{
		printf("Missing input mesh .obj \n");
		printf("Usage: %s [.obj file] [output_obj name] [max_impulse_value] [step]. \n", argv[0]);
		return 1;
	}

	printf("Starting application. \n");
	configFilesDir = string("configFiles\\");
	configFilename = string(argv[1]);
	outputFilename = string(argv[2]);
	sscanf(argv[3], "%lf", &impulse);
	sscanf(argv[4], "%lf", &max_impulse);
	sscanf(argv[5], "%lf", &step);

	if (argc == 7) //If there is a tool mesh obj
		extraSceneGeometryFilename = string(argv[6]);
	else
		extraSceneGeometryFilename = "__none";

	printf("Loading file %s \n", configFilename.c_str());

	initConfigurations();
	

	initGLUT(argc, argv, windowTitleBase, windowWidth, windowHeight, &windowID);
	initGraphics(windowWidth, windowHeight); // more OpenGL initialization calls

	initGLUI();
	initScene();

	//Start UI Loop
	glutMainLoop();

	return 0;
}