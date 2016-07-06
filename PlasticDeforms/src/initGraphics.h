#ifndef _INITGRAPHICS_H_
#define _INITGRAPHICS_H_

#ifdef WIN32
	#include <windows.h>
#endif

#include "openGL-headers.h"
#include "GL/glui.h"
#include "camera.h"

void initGLUT(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID);

void initCamera(double cameraRadius,
	double cameraLongitude, double cameraLattitude,
	double focusPosX, double focusPosY, double focusPosZ,
	double camera2WorldScalingFactor,
	double * zNear, double * zFar,
	SphericalCamera ** camera);

void initGraphics(int windowWidth, int windowHeight);

#endif