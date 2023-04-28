#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <chrono>
#include <limits>
#include <cmath>

#include "GLHeaders.h"

#define DPI_SCALE 1.0f
#define MSAA_SAMPLES 16

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "ShaderProgram.h"
#include "Buffer.h"
#include "GridSolver.h"
#include "Vertex.h"
#include "RecordVideo.hpp"
#include "ArcBall.h"

int winWidth = 1280, winHeight = 960;
bool simulating = false;
int fps = 30; // dipslay FPS (not simulation)

ImageRecorder imageRecorder;
int speed = 4;
bool recordFrames = false; // toggled by keyboard
bool recordFrame = false; // set to record when stepped in display`
//bool showFilament = false;
bool showGrid = false;
bool showVel = true;
bool showTemp = true;

bool velMode_all = true, velMode_long = false, velMode_lat = false, velMode_ver = false;
int velIdx_long = 0, velIdx_lat = 0, velIdx_layer = 0;
bool tpMode_all = true, tpMode_mer = false, tpMode_zon = false, tpMode_ver = false;
int tpIdx_long = 0, tpIdx_lat = 0, tpIdx_layer = 0;

//float velFactor = 1000.f;
//float srcTemp = 50.f;
//float pointedTemp = std::numeric_limits<double>::quiet_NaN();

// camera controls
/* arc ball */
ArcBall arcBall;
/* button state */
int buttonState = 0;
/* distance to object */
float zoom = 1;
/* distance to object min */
const float minDistance = .25;
/* distance to object max */
const float maxDistance = 16;

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::shared_ptr;
using glm::mat4;
using glm::vec3;
using std::chrono::steady_clock;

void keyboardCB(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS) {
		switch (key) {
            case GLFW_KEY_ESCAPE:
                glfwSetWindowShouldClose(window, GLFW_TRUE);
                break;
            case GLFW_KEY_S:
                step();
                //updateGrid();
                if (recordFrames)
                    recordFrame = true;
                break;
            case GLFW_KEY_R:
                initGrid();
                break;
            case GLFW_KEY_V:
                recordFrames = !recordFrames;
                break;
            case GLFW_KEY_SPACE:
                simulating = !simulating;
                break;
            case GLFW_KEY_UP:
                if (buttonState == 0)
                    arcBall.R = (mods == GLFW_MOD_SHIFT) ?
                                glm::rotate(glm::mat4(1.0f), (float) -(.5f * M_PI), glm::vec3(1, 0, 0)) :
                                glm::rotate(glm::mat4(1.f), (float) (.5f * M_PI), glm::vec3(1, 0, 0)) * arcBall.R;
                break;
            case GLFW_KEY_LEFT:
                if (buttonState == 0) {
                    if (mods == GLFW_MOD_SHIFT)
                        arcBall.R = glm::rotate(glm::mat4(1.0f), (float) -(.5f * M_PI), glm::vec3(1, 0, 0));
                    arcBall.R = glm::rotate(glm::mat4(1.f), (float) (.5f * M_PI), glm::vec3(0, 1, 0)) * arcBall.R;
                }
                break;
            case GLFW_KEY_RIGHT:
                if (buttonState == 0) {
                    if (mods == GLFW_MOD_SHIFT)
                        arcBall.R = glm::rotate(glm::mat4(1.f), (float) -(.5f * M_PI), glm::vec3(1, 0, 0));
                    arcBall.R = glm::rotate(glm::mat4(1.f), (float) -(.5f * M_PI), glm::vec3(0, 1, 0)) * arcBall.R;
                }
                break;
            case GLFW_KEY_DOWN:
                if (buttonState == 0)
                    arcBall.R = (mods == GLFW_MOD_SHIFT) ?
                                glm::rotate(glm::mat4(1.f), (float) M_PI, glm::vec3(1, 0, 0)) :
                                glm::rotate(glm::mat4(1.f), (float) -(.5f * M_PI), glm::vec3(1, 0, 0)) * arcBall.R;
                break;
        }
    }
}

//void getIJ(double x, double y, int &i, int &j) {
//	y = winHeight - y;
//	double x_offset, y_offset;
//	x_offset = (winWidth - frameWidth) / 2;
//	y_offset = (winHeight - frameHeight) / 2;
//	j = (int)((x - x_offset) / cellSize) + 1;
//	i = (int)((y - y_offset) / cellSize) + 1;
//}

void cursorPosCB(GLFWwindow* window, double xpos, double ypos) {
//	static double xprv = 0, yprv = 0;
	if (!ImGui::GetIO().WantCaptureMouse) {
//		int i, j;
//		getIJ(xpos, ypos, i, j);
//		if (0 < i && i <= longs && 0 < j && j <= lats)
//			pointedTemp = cur_tp[i][j];
//		else
//			pointedTemp = std::numeric_limits<double>::quiet_NaN();
//		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
//			double dx = xpos - xprv;
//                        double dy = ypos - yprv;
//			if (0 < i && i <= longs && 0 < j && j <= lats) {
//				double vx = velFactor * +dx;
//				double vy = velFactor * -dy;
//				vxSrc[i][j] = vx;
//				vySrc[i][j] = vy;
//			}
//		}
        if (buttonState == 1) {
            double x;
            double y;
            int windowWidth, windowHeight;
            glfwGetWindowSize(window, &windowWidth, &windowHeight);
            glfwGetCursorPos(window, &x, &y);
            arcBall.updateRotation(x, y, windowWidth, windowHeight);
        }
	}
//	xprv = xpos;
//	yprv = ypos;
}

void mouseButtonCB(GLFWwindow* window, int button, int action, int mods)
{
	if (!ImGui::GetIO().WantCaptureMouse) {
		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_RELEASE) {
                buttonState = 0;
                return;
            }
            double x;
            double y;
            int windowWidth, windowHeight;
            glfwGetWindowSize(window, &windowWidth, &windowHeight);
            glfwGetCursorPos(window, &x, &y);
            arcBall.startRotation(x, y, windowWidth, windowHeight);
            buttonState = 1;
		}
//        else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
//			if (action == GLFW_PRESS) {
//				double x, y;
//				glfwGetCursorPos(window, &x, &y);
//				x = x - (winWidth-frameWidth)/2;
//				y = (frameHeight+winHeight)/2 - y;
//				if (0 < x && x < frameWidth && 0 < y && y < frameHeight)
//					addFilament(x, y);
//			}
//		}
	}
}

void scrollCB(GLFWwindow* window, double xoffset, double yoffset) {
	if (!ImGui::GetIO().WantCaptureMouse) {
        zoom = glm::clamp((float)(zoom * ::pow(2, yoffset)), minDistance, maxDistance);

//        double x, y;
//		int i, j;
//		glfwGetCursorPos(window, &x, &y);
//		getIJ(x, y, i, j);
//		if (0 < i && i <= rows && 0 < j && j <= cols) {
//			tpSrc[i][j] += 0.1 * yoffset;
//			pointedTemp = tpSrc[i][j];
//		}
	}
}

int main(void) {
	GLFWwindow* window;

	if (!glfwInit()) {
		cerr << "glfwInit() failed" << endl;
		exit(EXIT_FAILURE);
	}

	// GL 4.0
	// GLSL 400
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

	glfwWindowHint(GLFW_SAMPLES, MSAA_SAMPLES);
	window = glfwCreateWindow(winWidth, winHeight, "Fluid", NULL, NULL);
	if (!window) {
		cerr << "glfwCreateWindow() failed" << endl;
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwSetKeyCallback(window, keyboardCB);
	glfwSetCursorPosCallback(window, cursorPosCB);
	glfwSetMouseButtonCallback(window, mouseButtonCB);
	glfwSetScrollCallback(window, scrollCB);
	glfwMakeContextCurrent(window);

	GLenum errn = glewInit();
	if (GLEW_OK != errn) {
		cerr << "glewInit() failed: " << glewGetErrorString(errn) << endl;
		exit(EXIT_FAILURE);
	}

	glEnable(GL_MULTISAMPLE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	ShaderProgram shader("../resources/shaders/simple.vs", "../resources/shaders/simple.fs");

	vector<size_t> attrSizes{ 3, 4 };

	initGrid();

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsLight();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 400");
#ifdef DPI_SCALE
	ImFontConfig font_cfg;
	font_cfg.SizePixels = (int)(13 * DPI_SCALE);
	ImGui::GetIO().Fonts->AddFontDefault(&font_cfg);
	ImGui::GetStyle().ScaleAllSizes(DPI_SCALE);
#endif

	Vertex frameVertices[2 * (longs + lats + 2) * (layers + 1)];
	size_t k = 0;
	const float c = 0.2f;
    for (size_t l = 0; l <= layers; ++l) {
        for (size_t i = 0; i <= longs; ++i) {
            frameVertices[k++] = Vertex((float)i * cellSize, (float)0, (float)l * cellSize, 0.2f, 1.0f, 0.2f, c);
            frameVertices[k++] = Vertex((float)i * cellSize, (float)frameDepth, (float)l * cellSize, 1.0f, 1.0f, 1.0f, c);
        }
        for (size_t j = 0; j <= lats; ++j) {
            frameVertices[k++] = Vertex((float)0, (float)j * cellSize, (float)l * cellSize, 0.2f, 1.0f, 0.2f, c);
            frameVertices[k++] = Vertex((float)frameWidth, (float)j * cellSize, (float)l * cellSize, 1.0f, 1.0f, 1.0f, c);
        }
    }
	Buffer frameBuffer(attrSizes, (GLfloat*)frameVertices, sizeof(frameVertices));

	auto t_now = steady_clock::now();
	while (!glfwWindowShouldClose(window)) {

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::Begin("Sim Setting");
		ImGui::Text("Shortcut:\n\tspace - run/stop\n\ts     - step\n\tr     - restart\n\tv     - record\n\tESC   - exit");
		ImGui::Text("\tright button drag\t- rotate around\n\tscroll\t- zoom\n\tarrow keys\t- 90 degree rotate\n\tshift + arrow keys\t- preset rotations");
		ImGui::Text("simTime: %6.2f", simTime);
//		ImGui::Text("srcTemp: %6.2f", pointedTemp);
		ImGui::Text("avgTemp: %6.2f", (float)getReferenceTemperature(*cur_tp));
		ImGui::SliderInt("speed", &speed, 1, 10);
		ImGui::Checkbox("recordFrames", &recordFrames);
		ImGui::Checkbox("useCoriolisForce", &useCoriolis);
		ImGui::SliderFloat("dt", &dt, 0.0001f, 0.1f);
		ImGui::SliderFloat("dx", &dx, 0.01f, 5.0f);
		ImGui::SliderFloat("nv", &nv, 0.f, 20.f);
		ImGui::SliderFloat("kappa", &kappa, 0.f, 5.f);
		ImGui::SliderFloat("buoyancy", &buoyancy, 0.f, 5.f);
		ImGui::SliderFloat("topTemp", &top_temp, -10.f, 10.f);
		ImGui::SliderFloat("botTemp", &bot_temp, -10.f, 10.f);
		ImGui::SliderFloat("planetaryRotation", &planetary_rotation, 0.f, 10.f);
		ImGui::SliderFloat("externalHeatFactor", &external_heat_factor, 0.f, 1.0f);
		ImGui::SliderInt("iterations", &iterations, 1, 100);
		ImGui::End();

		ImGui::Begin("Viz Setting");
//		ImGui::Checkbox("showFilament", &showFilament);
        ImGui::Checkbox("showGrid", &showGrid);
        ImGui::Checkbox("showVel", &showVel);
        ImGui::Checkbox("showTemp", &showTemp);
        if (ImGui::TreeNode("Velocity")){
            ImGui::SliderFloat("velScale", &velScale, 1.f, 200.f);
            ImGui::Checkbox("All", &velMode_all);
            if (velMode_all) ImGui::BeginDisabled();
            ImGui::Checkbox("Longitude", &velMode_long);
            ImGui::SliderInt("X Index", &velIdx_long, 0, longs - 1);
            ImGui::Checkbox("Latitude", &velMode_lat);
            ImGui::SliderInt("Y Index", &velIdx_lat, 0, lats - 1);
            ImGui::Checkbox("Height", &velMode_ver);
            ImGui::SliderInt("Z Index", &velIdx_layer, 0, layers - 1);
            if (velMode_all) ImGui::EndDisabled();
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Temperature")){
            ImGui::Checkbox("All", &tpMode_all);
            if (tpMode_all) ImGui::BeginDisabled();
            ImGui::Checkbox("Longitude", &tpMode_mer);
            ImGui::SliderInt("X Index", &tpIdx_long, 0, longs);
            ImGui::Checkbox("Latitude", &tpMode_zon);
            ImGui::SliderInt("Y Index", &tpIdx_lat, 0, lats);
            ImGui::Checkbox("Height", &tpMode_ver);
            ImGui::SliderInt("Z Index", &tpIdx_layer, 0, layers);
            if (tpMode_all) ImGui::EndDisabled();
            ImGui::TreePop();
        }
        ImGui::End();

		ImGui::Render();

		glfwGetFramebufferSize(window, &winWidth, &winHeight);

		glViewport(0, 0, winWidth, winHeight);
		glClear(GL_COLOR_BUFFER_BIT);

		mat4 pv = glm::ortho(-(float)winWidth / 2, (float)winWidth / 2, -(float)winHeight / 2, (float)winHeight / 2, -(float)4*frameWidth, (float)4*frameWidth); // Note that glm::ortho is a template function. Must cast the arguments to float to use the function with the correct type!
        mat4 m = glm::translate(vec3(-frameWidth/2, -frameDepth/2, -frameHeight / 2));
//		mat4 m_zoom = glm::translate(vec3(0.f, 0.f, -zoom * cellSize));
		mat4 m_zoom = glm::scale(vec3(zoom));
		shader.use();
		shader.set_mat4("MVP", pv * m_zoom * arcBall.R * m);
//		shader.set_mat4("MVP", pv * arcBall.R * m);

		if (simulating && std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now() - t_now).count() > 1000 / fps) {
			for (int i = 0; i < speed; ++i)
				step();
			updateGrid();
			t_now = steady_clock::now();
			if (recordFrames)
				recordFrame = true;
		}
		//updateGrid();

//		Buffer gridBuffer(attrSizes, (GLfloat*)gridVertices, sizeof(gridVertices), gridIndices, sizeof(gridIndices));
        if (showTemp) {
            Buffer gridBuffer(attrSizes, (GLfloat*) &gridVertices.buffer[0], (longs + 1) * (lats + 1) * (layers + 1) * sizeof (Vertex));
            if (tpMode_all){
                glPointSize(5.f);
                gridBuffer.draw(GL_POINTS);
            } else {
                gridBuffer.bind();
                if (tpMode_mer)
                    glDrawElements(GL_TRIANGLES, 6 * lats * layers, GL_UNSIGNED_INT, &gridIndices_long[tpIdx_long][0]);
                if (tpMode_zon)
                    glDrawElements(GL_TRIANGLES, 6 * longs * layers, GL_UNSIGNED_INT, &gridIndices_lat[tpIdx_lat][0]);
                if (tpMode_ver)
                    glDrawElements(GL_TRIANGLES, 6 * longs * lats, GL_UNSIGNED_INT, &gridIndices_vert[tpIdx_layer][0]);
                gridBuffer.unbind();
            }
        }


        if (showVel) {
            Buffer velBuffer(attrSizes, (GLfloat*) &velVertices.buffer[0][0], 2 * longs * lats * layers * sizeof (Vertex));
            if (velMode_all) {
                velBuffer.draw(GL_LINES);
            } else {
                velBuffer.bind();
                if (velMode_long)
                    glDrawElements(GL_LINES, 2 * lats * layers, GL_UNSIGNED_INT, &velIndices_long[velIdx_long][0]);
                if (velMode_lat)
                    glDrawElements(GL_LINES, 2 * longs * layers, GL_UNSIGNED_INT, &velIndices_lat[velIdx_lat][0]);
                if (velMode_ver)
                    glDrawElements(GL_LINES, 2 * longs * lats, GL_UNSIGNED_INT, &velIndices_vert[velIdx_layer][0]);
                velBuffer.unbind();
            }
        }

		if (showGrid)
			frameBuffer.draw(GL_LINES);

		if (recordFrame)
		{
			recordFrame = false;
			imageRecorder.writeCurrentFrameToFile(window);
		}

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);

	glfwTerminate();
	exit(EXIT_SUCCESS);
}
