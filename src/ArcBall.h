// Jonathan McIntosh 260850688

#pragma once
#ifndef ARCBALL_H
#define ARCBALL_H

#include <string>
#include <vector>

#include "GLHeaders.h"

/**
 * By Alexandre Mercier-Aubin
 * This class manages the camera interaction and motions using an ArcBall.
 */
class ArcBall
{
public:
	ArcBall();
	virtual ~ArcBall();

	// this starts the tracking motion according to inital mouse positions. It should be called on mouse down.
	void startRotation(double mousex, double mousey, int windowWidth, int windowHeight);

	// This generates a rotation based on the mouse displacement on screen with respect to the initial mouse position on click. Call this on mouse move.
	void updateRotation(double mousex, double mousey, int windowWidth, int windowHeight);

	// This dun
	glm::vec3 computeVecFromMousePos(double mousex, double mousey, int windowWidth, int windowHeight);

	// the arcball generated rotation matrix to change the camera's view
	glm::mat4 R;

	// parameters to tune the speed of rotation
	double fit;
	double gain;

	const double eps = 1e-2;

private:
	// Stores the initial rotation matrix, prior to the update
	glm::mat4 Rmem;

	// stores the pre initial projected vec
	glm::vec3 p0;

	// stores the current projected vec
	glm::vec3 p1;
};

#endif
