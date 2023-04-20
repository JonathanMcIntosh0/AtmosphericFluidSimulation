// Jonathan McIntosh 260850688

#include "ArcBall.h"
//#include "MatrixStack.h"
//#include "Program.h"

#include <iostream>
#include <cassert>

//#include "GLSL.h"

using namespace std;

ArcBall::ArcBall() : R(glm::mat4(1.0)), Rmem(glm::mat4(1.0)), p0(glm::vec3(1.0)), p1(glm::vec3(1.0)), fit(0.5), gain(5.0)
{
}

ArcBall::~ArcBall()
{
}

glm::vec3 ArcBall::computeVecFromMousePos(double mousex, double mousey, int windowWidth, int windowHeight)
{
	double smallestDimension = (windowWidth > windowHeight ? windowHeight : windowWidth);
	double radius = smallestDimension / fit;

	glm::vec3 v = glm::vec3(0, 0, 0);
	v.x = (mousex - 0.5 * windowWidth) / radius;
	v.y = -(mousey - 0.5 * windowHeight) / radius;

	double r = v.x * v.x + v.y * v.y;
	if (r > 1.0)
	{
		double s = 1.0 / sqrt(r);
		v.x = s * v.x;
		v.y = s * v.y;
		v.z = 0;
	}
	else
	{
		v.z = sqrt(1.0 - r);
	}
	return v;
}

double computeVectorAngle(glm::vec3 &v1, glm::vec3 &v2)
{
	double vDot = glm::dot(v1, v2);
	if (vDot < -1.0)
		vDot = -1.0;
	if (vDot > 1.0)
		vDot = 1.0;
	return ((double)(acos(vDot)));
}

void ArcBall::startRotation(double mousex, double mousey, int windowWidth, int windowHeight)
{
	Rmem = R;
	p0 = computeVecFromMousePos(mousex, mousey, windowWidth, windowHeight);
}

void ArcBall::updateRotation(double mousex, double mousey, int windowWidth, int windowHeight)
{
	R = Rmem;
	glm::vec3 p1 = computeVecFromMousePos(mousex, mousey, windowWidth, windowHeight);
	glm::vec3 pCross = glm::cross(p0, p1);
	float angle = computeVectorAngle(p0, p1) * gain;
	if (angle > eps) {
		glm::mat4 rotationUpdate = glm::rotate(glm::mat4(1.0), angle, pCross);
		R = rotationUpdate * R;
	}
}
