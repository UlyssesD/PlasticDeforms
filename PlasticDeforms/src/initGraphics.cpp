#ifndef M_PI
	#define M_PI 3.141592653589793238462643
#endif

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include "initGraphics.h"

extern void displayFunction(void);
extern void keyboardFunction(unsigned char key, int x, int y);
extern void mouseButtonActivityFunction(int button, int state, int x, int y);
extern void mouseMotionFunction(int x, int y);
extern void idleFunction(void);
extern void reshape(int, int);

// initialize GLUT
void initGLUT(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID)
{
	// Initialize GLUT.
	glutInit(&argc, argv);
	//glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(windowWidth, windowHeight);
	*windowID = glutCreateWindow(windowTitle);

	// Setup GLUT callbacks.
	glutDisplayFunc(displayFunction);
	glutMouseFunc(mouseButtonActivityFunction);
	glutKeyboardFunc(keyboardFunction);
	//GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
	//GLUI_Master.set_glutSpecialFunc(specialFunction);
	//GLUI_Master.set_glutReshapeFunc(reshape);

	glutMotionFunc(mouseMotionFunction);
	glutIdleFunc(idleFunction);
}

void initCamera(double cameraRadius, double cameraLongitude, double cameraLattitude, double focusPosX, double focusPosY, double focusPosZ, double camera2WorldScalingFactor, double * zNear, double * zFar, SphericalCamera ** camera)
{
	double focusPos[3] = { focusPosX, focusPosY, focusPosZ };

	*zNear = cameraRadius * 0.01;
	*zFar = cameraRadius * 100;

	double upPos[3] = { 0,1,0 };
	*camera = new SphericalCamera(cameraRadius, 1.0 * cameraLongitude / 360 * (2 * M_PI), 1.0 * cameraLattitude / 360 * (2 * M_PI), focusPos, upPos, 0.05, camera2WorldScalingFactor);

	(*camera)->SetOrigin(focusPos);
}

void initGraphics(int windowWidth, int windowHeight)
{
	// clear to white
	//glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

	// clear to light blue
	//glClearColor(233.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

	// clear to gray
	//glClearColor(196.0 / 256, 196.0 / 256, 196.0 / 256, 0.0);

	// clear to brown
	//glClearColor(255.0 / 256, 156.0 / 256, 17.0 / 256, 0.0);

	// clear to medical blue
	glClearColor(148.0 / 256, 199.0 / 256, 211.0 / 256, 0.0);

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_STENCIL_TEST);
	//glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	reshape(windowWidth, windowHeight);

	printf("Graphics initialization complete.\n");
}