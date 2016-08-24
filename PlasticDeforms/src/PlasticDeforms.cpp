#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <limits>
#include <direct.h>
#include <vector>
#include <map>
#include <list>

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

//collisionDetection inclusions
#include "CTCD.h"
#include <Eigen/Core>

//Project inclusions
#include "initGraphics.h"

using namespace std;

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
double threshold = 1E-5;
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

//map<int, Eigen::Vector3d> vtxs;
map<int, pair<Eigen::Vector3d, double>> vtxs;
//map<int, Vec3d> vtxs;

//Vec3d trajectory[4] = { Vec3d(-0.05, 0.03, 0), Vec3d(-0.05, 0.015, 0), Vec3d(0.05, 0.015, 0), Vec3d(0.05, 0.03, 0) };
//Vec3d trajectory[3] = { Vec3d(-0.05, 0.03, 0), Vec3d(-0.05, 0.015, 0), Vec3d(0.05, 0.015, 0) };
//Vec3d trajectory[10] = { Vec3d(-0.08, 0.03, 0), Vec3d(-0.08, 0.015, 0), Vec3d(-0.06, 0.015, 0.02), Vec3d(-0.04, 0.015, 0.04),  Vec3d(-0.02, 0.015, 0.06), Vec3d(0.0, 0.015, 0.08), Vec3d(0.02, 0.015, 0.06), Vec3d(0.04, 0.015, 0.04), Vec3d(0.06, 0.015, 0.02), Vec3d(0.08, 0.015, 0.0) };
//for point_tool_x
Vec3d trajectory[10] = { Vec3d(-0.08, 0.03, 0), Vec3d(-0.08, 0.01, 0), Vec3d(-0.06, 0.01, 0.02), Vec3d(-0.04, 0.01, 0.04),  Vec3d(-0.02, 0.01, 0.06), Vec3d(0.0, 0.01, 0.08), Vec3d(0.02, 0.01, 0.06), Vec3d(0.04, 0.01, 0.04), Vec3d(0.06, 0.01, 0.02), Vec3d(0.08, 0.01, 0.0) };
//Vec3d trajectory[3] = { Vec3d(-0.08, 0.03, 0), Vec3d(-0.08, 0.012, 0), Vec3d(0.08, 0.012, 0) };

//list of faces under analysis
list<int> faces;
void applyImpulseForce()
{
	cout << "Applying an impulse force of" << impulse << "Newton along the Y axes on selected vertices" << endl;
	
	double externalForce[3] = { -1*impulse/2, -1*impulse/2, 0 };

	//printf("Externalforce: [%f, %f, %f]\n", externalForce[0], externalForce[1], externalForce[2]);

	double cos = acos(dot(Vec3d(0, impulse, 0), Vec3d(0, 1, 0))/ len(Vec3d(0, impulse, 0)) * len(Vec3d(0, 1, 0)));
	//printf("Is force perpendicular: %d. \n", cos == PI);

	for (int i = 0; i < numverts; i++) {
		//printf("Applying an impulse force of %f Newton along the Y axes on vertex %d. \n", impulse, vertices[i]);
		//	renderingModalMatrix->ProjectSingleVertex(vertices[i], externalForce[0], externalForce[1], externalForce[2], fq);
			 
		//	for (int j = 0; j < r; j++)
		//		fqBase[j] = fqBase[j] + fq[j];
	
	renderingModalMatrix->ProjectSingleVertex(vertex,
		externalForce[0], externalForce[1], externalForce[2], fq);
	}

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
				cout << "The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep" << endl;
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

	cout << "System kinetic energy: " << implicitNewmarkDense->GetKineticEnergy() << endl;
	
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

vector<double> getPositionVector(SceneObjectDeformable* obj) 
{
	vector<double> pos;

	for (int j = 0; j < obj->GetNumVertices(); j++)
	{

		double p[3];
		obj->GetSingleVertexPositionFromBuffer(j, &p[0], &p[1], &p[2]);
		for (int i = 0; i<3; i++)
		{
			pos.push_back(p[i]);
		}
	}

	return pos;
}

vector<int> getFacesVector(SceneObjectDeformable* obj) 
{
	vector<int> faceidx;

	for (int j = 0; j < obj->GetNumFaces(); j++)
	{

		ObjMesh::Face t = obj->GetMesh()->getGroup("default").getFace(j);
		for (int i = 0; i<3; i++)
		{

			ObjMesh::Vertex v = t.getVertex(i);
			faceidx.push_back(v.getPositionIndex());
		}
	}

	return faceidx;
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

void simulateImpact() 
{
	//double externalForce[3] = { -1*impulse, -1 * impulse, 0 };

	//printf("Externalforce: [%f, %f, %f]\n", externalForce[0], externalForce[1], externalForce[2]);

	//double cos = acos(dot(Vec3d(0, impulse, 0), Vec3d(0, 1, 0)) / len(Vec3d(0, impulse, 0)) * len(Vec3d(0, 1, 0)));
	//printf("Is force perpendicular: %d. \n", cos == PI);

	//for (std::map<int, Vec3d>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
	//for (std::map<int,Eigen::Vector3d>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
	for (std::map<int, pair<Eigen::Vector3d, double>>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
	{
		//double externalForce[3] = { -1 * impulse * it->second[0], -1 * impulse * it->second[1], -1 * impulse * it->second[2] };
		double externalForce[3] = { -1 * impulse * it->second.first[0], -1 * impulse * it->second.first[1], -1 * impulse * it->second.first[2] };
		renderingModalMatrix->ProjectSingleVertex(it->first, externalForce[0], externalForce[1], externalForce[2], fq);

		for (int j = 0; j < r; j++)
			fqBase[j] = fqBase[j] + fq[j];
	}
	
	vtxs.clear();
/*
	for (int i = 0; i < numverts; i++) {
		//printf("Applying an impulse force of %f Newton along the Y axes on vertex %d. \n", impulse, vertices[i]);
		renderingModalMatrix->ProjectSingleVertex(vertices[i], externalForce[0], externalForce[1], externalForce[2], fq);

		for (int j = 0; j < r; j++)
			fqBase[j] = fqBase[j] + fq[j];
	}
*/
	//renderingModalMatrix->ProjectSingleVertex(vertex,
	//	externalForce[0], externalForce[1], externalForce[2], fq);

	memcpy(fq, fqBase, sizeof(double) * r);
	// set the reduced external forces
	implicitNewmarkDense->SetExternalForces(fq);



	//printf("System kinetic energy: %f \n", implicitNewmarkDense->GetKineticEnergy());

//	do {

		for (int i = 0; i < substepsPerTimeStep; i++)
		{
			int code = implicitNewmarkDense->DoTimestep();
			if (code != 0)
			{
				//printf("The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep.\n");
				implicitNewmarkDense->ResetToRest();
			}

			//memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

			// compute u=Uq
			deformableObjectRenderingMeshReduced->Setq(implicitNewmarkDense->Getq());
			deformableObjectRenderingMeshReduced->Compute_uUq();

			deformableObjectRenderingMeshReduced->BuildNormals();

			displayFunction();

		}

//	} while (implicitNewmarkDense->GetKineticEnergy() > threshold);
/*
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
	*/
}

void testCollisions(Eigen::Vector3d dist, Vec3d d) 
{
	//compute center of mass of tool
	//Vec3d c = extraSceneGeometry->GetMesh()->computeCenterOfMass_Vertices();
	//Eigen::Vector3d com(c[0], c[1], c[2]);

	//numverts = 0;
	//vertices = (int *)malloc(deformableObjectRenderingMeshReduced->GetNumVertices() * sizeof(int));
	ObjMesh* surface_mesh = deformableObjectRenderingMeshReduced->GetMesh();
	//surface_mesh->initTriangleLookup();

	//ObjMesh* tool_mesh = extraSceneGeometry->GetMesh();
	//tool_mesh->initTriangleLookup();

	//compute collisions by checking surface vertex - tool triangles collisions
	//for(int i = 0; i < surface_mesh->getNumFaces(); i++)
	//for (int i = 0; i < extraSceneGeometry->GetNumFaces(); i++)
	for(list<int>::iterator itFaces = faces.begin(); itFaces != faces.end(); ++ itFaces)
	{
		/*
		//get i-th triangle 
		int triangle[3];
		tool_mesh->getTriangle(i, &triangle[0], &triangle[1], &triangle[2]);
		*/
		int triangle[3];
		//surface_mesh->getTriangle(i, &triangle[0], &triangle[1], &triangle[2]);
		surface_mesh->getTriangle(*itFaces, &triangle[0], &triangle[1], &triangle[2]);

		//get position values for i-th triangle's vertices
		double pos[9];
		for (int t = 0; t < 3; t++)
			deformableObjectRenderingMeshReduced->GetSingleVertexPosition(triangle[t], &pos[3 * t + 0], &pos[3 * t + 1], &pos[3 * t + 2]);
			//extraSceneGeometry->GetSingleVertexPositionFromBuffer(triangle[t], &pos[3 * t + 0], &pos[3 * t + 1], &pos[3 * t + 2]);

		Eigen::Vector3d t0(pos[0], pos[1], pos[2]);
		Eigen::Vector3d t1(pos[3], pos[4], pos[5]);
		Eigen::Vector3d t2(pos[6], pos[7], pos[8]);

		for (int j = 0; j< extraSceneGeometry->GetNumVertices(); j++ )
		//for(int j = 0; j < deformableObjectRenderingMeshReduced->GetNumVertices(); j++)
		{
			double vert[3];
			//deformableObjectRenderingMeshReduced->GetSingleVertexPosition(j, &vert[0], &vert[1], &vert[2]);
			extraSceneGeometry->GetSingleVertexPositionFromBuffer(j, &vert[0], &vert[1], &vert[2]);
			Eigen::Vector3d v(vert[0], vert[1], vert[2]);

			double eta = 1e-4;
			double t = 0;
			
			//bool collision = CTCD::vertexFaceCTCD(v, t0 - dist, t1 - dist, t2 - dist, v, t0, t1, t2, eta, t);
			bool collision = CTCD::vertexFaceCTCD(v - dist, t0, t1, t2, v, t0, t1, t2, eta, t);

			if (collision)
			{
				//Calculate (normalized) face centroid and force direction
				
				//Vec3d centroid = Vec3d(pos[0], pos[1], pos[2]) + Vec3d(pos[3], pos[4], pos[5]) + Vec3d(pos[6], pos[7], pos[8]) / 3;
				//Vec3d force = (Vec3d(vert[0], vert[1], vert[2]) - d) - centroid;
				
				Eigen::Vector3d force = (v - dist) - ((t0 + t1 + t2) / 3);
				//Eigen::Vector3d force =  (v) - ((t0 + t1 + t2) / 3);

				//force direction based on center of mass
				//Eigen::Vector3d force = (com) - ((t0 + t1 + t2) / 3);

				//printf("Collision between vertex %d and face [%d, %d, %d]!\n", j, triangle[0], triangle[1], triangle[2]);
				double l = force.norm();
				force.normalize();

				for (int t = 0; t < 3; t++)
				{
					map<int, pair<Eigen::Vector3d, double>>::iterator it = vtxs.find(triangle[t]);
					if (it == vtxs.end())
						vtxs.emplace(triangle[t], pair<Eigen::Vector3d, double>(force, l));
					else if(l < it->second.second)
					{
						it->second.first = force;
						it->second.second = l;
					}
					//vertices[numverts] = triangle[t];
					//numverts++;
				}

				/*
				vertices[numverts] = j;
				numverts++;
				*/
			}
		}
	}
	if (vtxs.size() != 0)
	{
		cout << "Applying forces on vertices: ";
		//for (std::map<int, Vec3d>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
		//for (std::map<int,Eigen::Vector3d>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
		for (std::map<int, pair<Eigen::Vector3d, double>>::iterator it = vtxs.begin(); it != vtxs.end(); ++it)
			cout <<  it->first << " ";
		cout << endl;

	}
	//free(vertices);
}

void findTopSurfaceVertices() {
	ObjMesh* surface_mesh = deformableObjectRenderingMeshReduced->GetMesh();

	for (unsigned int i = 0; i < surface_mesh->getNumGroups(); i++)
	{
		for (unsigned int iFace = 0; iFace < surface_mesh->getGroupHandle(i)->getNumFaces(); iFace++)
		{

			ObjMesh::Face * faceHandle = (ObjMesh::Face *) surface_mesh->getGroupHandle(i)->getFaceHandle(iFace); // get face whose number is iFace

			if (faceHandle->getNumVertices() < 3)
				cout << "Warning: encountered a face (group = " << i << ", face = " << iFace << ") with fewer than 3 vertices." << endl;

			Vec3d normal = surface_mesh->computeFaceNormal(*faceHandle);

			if (normal[1] > 0 && normal[0] == 0 && normal[2] == 0)
				faces.push_front(iFace);
		}

	}
	/*
	cout << "Number of total faces: " << surface_mesh->getNumFaces() << "; \t size of faces list: " << faces.size() << "; \t top surface faces: ";
	for (list<int>::iterator it = faces.begin(); it != faces.end(); ++it)
		cout << *it << " ";
	cout << endl;
	*/
}

void keyboardFunction(unsigned char key, int x, int y)
{
	Vec3d u; Vec3d  dist; Vec3d tool;
	double span;

	int subs = 20;

	//vector<int> mesh_faces; vector<int> tool_faces;
	//vector<double> mesh_pos; vector<double> tool_pos;

	switch (key)
	{
		case 27:
		case 'q':
			exit(0);
			break;
		case 'b':
			findTopSurfaceVertices();
			break;
		case 'f':
			cout << "Simulate tool swipe along plate." << endl;
			
			//distance of swipe
			dist = (trajectory[2] - trajectory[1]);
			span = dist[0] / 50;
			tool = trajectory[1];
			for (int substeps = 0; substeps < 50; substeps++)
			{
				tool[0] = tool[0] + span;
				tool.print();
				numverts = 0;
				vertices = (int *)malloc(deformableObjectRenderingMeshReduced->GetNumVertices() * sizeof(int));
				double v[3];
				for (int i = 0; i < deformableObjectRenderingMeshReduced->GetNumVertices(); i++) {
					deformableObjectRenderingMeshReduced->GetSingleVertexPosition(i, &v[0], &v[1], &v[2]);
					if (v[0] >= tool[0] - 0.002 && v[0] <= tool[0] + 0.002  && v[2] >= tool[2] - 0.001 && v[2] <= tool[2] + 0.001) {
						//deformableObjectRenderingMeshReduced->HighlightVertex(i);
						//printf("Adding v[%d]: %f, %f, %f; \n", i, v[0], v[1], v[2]);
						vertices[numverts] = i;
						numverts++;
						//renderingModalMatrix->ProjectSingleVertex(i, externalForce[0], externalForce[1], externalForce[2], fq);
					}

				}
				cout << "final size of impact area: " << numverts << ";" <<endl;

				simulateImpact();
			}
			cout << "System kinetic energy: " << implicitNewmarkDense->GetKineticEnergy() << endl;
			while (implicitNewmarkDense->GetKineticEnergy() > threshold) {

				for (int i = 0; i < substepsPerTimeStep; i++)
				{
					int code = implicitNewmarkDense->DoTimestep();
					if (code != 0)
					{
						//printf("The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep.\n");
						implicitNewmarkDense->ResetToRest();
					}

					//memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

					// compute u=Uq
					deformableObjectRenderingMeshReduced->Setq(implicitNewmarkDense->Getq());
					deformableObjectRenderingMeshReduced->Compute_uUq();

					deformableObjectRenderingMeshReduced->BuildNormals();

					displayFunction();

				}
			}

			for (int i = 0; i < r; i++)
			{
				fq[i] = 0;
				fqBase[i] = 0;
			}
			break;

		case 'v':
			cout << "Calculating impact area;\n";
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
			cout << "Final size of impact area:" << numverts << ";\n";
			break;


		case ' ':
			applyImpulseForce();
			break;
		case 'p':
			deformableObjectRenderingMeshReduced->GetMesh()->triangulate();
			extraSceneGeometry->GetMesh()->triangulate();
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
			cout << "Starting tool movement simulation \n";
			
			findTopSurfaceVertices();
			//load face vectors for tool and simulation mesh
			//mesh_faces = getFacesVector(deformableObjectRenderingMeshReduced);
			//tool_faces = getFacesVector(extraSceneGeometry);

			cout << "Moving vertices. \n";
			for (int p = 0; p < (size(trajectory) - 1); p++) {
				dist = (trajectory[p + 1] - trajectory[p]); /// len(trajectory[1] - trajectory[0]);

				printf("d: [%f %f %f] \n", dist[0], dist[1], dist[2]);

				displacements = (double *)malloc(3 * extraSceneGeometry->GetNumVertices() * sizeof(double));
				for (int substeps = 0; substeps < subs; substeps++)
				{
					for (int i = 0; i < extraSceneGeometry->GetNumVertices(); i++)
					{
						displacements[i * 3 + 0] = dist[0] / subs;
						displacements[i * 3 + 1] = dist[1] / subs;
						displacements[i * 3 + 2] = dist[2] / subs;
					}
					extraSceneGeometry->AddVertexDeformations(displacements);
					displayFunction();

					testCollisions(Eigen::Vector3d(dist[0] / subs, dist[1] / subs, dist[2] / subs), Vec3d(dist[0] / subs, dist[1] / subs, dist[2] / subs));
					//mesh_pos = getPositionVector(deformableObjectRenderingMeshReduced);
					//tool_pos = getPositionVector(extraSceneGeometry);

					simulateImpact();	
				}

				free(displacements);
			}
			/*
			while (implicitNewmarkDense->GetKineticEnergy() > threshold) {

				for (int i = 0; i < substepsPerTimeStep; i++)
				{
					int code = implicitNewmarkDense->DoTimestep();
					if (code != 0)
					{
						//printf("The integrator went unstable. Reduce the timestep, or increase the number of substeps per timestep.\n");
						implicitNewmarkDense->ResetToRest();
					}

					//memcpy(q, implicitNewmarkDense->Getq(), sizeof(double) * r);

					// compute u=Uq
					deformableObjectRenderingMeshReduced->Setq(implicitNewmarkDense->Getq());
					deformableObjectRenderingMeshReduced->Compute_uUq();

					deformableObjectRenderingMeshReduced->BuildNormals();

					displayFunction();

				}

			}

			for (int i = 0; i < r; i++)
			{
				fq[i] = 0;
				fqBase[i] = 0;
			}
			*/
			faces.clear();
			printf("Execution Completed.\n");
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

	//configFile.addOptionOptional("plasticThreshold", &plasticThreshold, plasticThreshold);

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
	implicitNewmarkDense->SetPlasticThreshold(plasticThreshold);
	//implicitNewmarkDense->UsePlasticDeformations(1);
	//implicitNewmarkDense->SetInternalForceScalingFactor(frequencyScaling * frequencyScaling);
	deformableObjectRenderingMeshReduced->GetMesh()->triangulate();
	deformableObjectRenderingMeshReduced->GetMesh()->initTriangleLookup();
	deformableObjectRenderingMeshReduced->BuildNeighboringStructure();
	deformableObjectRenderingMeshReduced->BuildNormals(85.0);
	//implicitNewmarkDense->UseStaticSolver(1);
	free(massMatrix);

	if (strcmp(extraSceneGeometryFilename.c_str(),"__none") != 0)
	{
		extraSceneGeometry = new SceneObjectDeformable(extraSceneGeometryFilename.c_str());
		extraSceneGeometry->GetMesh()->triangulate();
		deformableObjectRenderingMeshReduced->GetMesh()->initTriangleLookup();
		extraSceneGeometry->BuildNeighboringStructure();

		extraSceneGeometry->BuildNormals(85.0);
	}
	else
		extraSceneGeometry = NULL;

	// compute lowest frequency of the system (smallest eigenvalue of K)
	double * Kq = (double*)malloc(sizeof(double) * r * r);
	double * zero = (double*)calloc(r, sizeof(double));
	stVKReducedStiffnessMatrix->Evaluate(zero, Kq);

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