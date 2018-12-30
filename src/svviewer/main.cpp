
//=============================================================
//  This file is part of the SymVox (Symmetry Voxelization) software
//  Copyright (C) 2016 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://vic.crs4.it
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//=============================================================


#include <fstream>
#include <string.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <symvox/encoded_octree.hpp>
#include <symvox/encoded_svdag.hpp>
#include <symvox/encoded_ussvdag.hpp>
#include <symvox/encoded_ssvdag.hpp>
#include "camera.hpp"
#include "renderer_monitor.hpp"
#include "octree_dda_renderer.hpp"

#define UPDATE_INFO_TIME 200.0 // ms
#define SCREEN_WIDTH	1280
#define SCREEN_HEIGHT	720

bool finish = false;
bool printCamera = false;
float frameTime = 0;

OctreeDDARenderer * renderer;
EncodedOctree* encoded_octree;
RendererMonitor::FrameStats frameStats;
GLFWwindow * window;
Camera * cam;

void printHelp() {
	printf(
		"* SymVox Viewer keys:\n"
		"\t- 1,2 : decrease / increase octree draw level\n"
		"\t- 3,4 : decrease / increase max traversal iterations\n"
		"\t- 5,6 : decrease / increase pixel tolerance\n"
		"\t- R   : switch render mode [ N.Steps | Depth | DrawLevel ]\n"
		"\t- O   : toggle beam optimization\n"
		"\t- ESC : exit\n"
		"\n");
	cam->printControls();
}

void optionsKeyCallback(GLFWwindow *win, int key, int scancode, int action, int mods) {
	if (action != GLFW_PRESS) return;

	// generic keys
	if (key == GLFW_KEY_ESCAPE) finish = true;
	if (key == GLFW_KEY_1) renderer->decDrawLevel();
	if (key == GLFW_KEY_2) renderer->incDrawLevel();
	if (key == GLFW_KEY_3) renderer->decGPUTraversalMaxIters(10);
	if (key == GLFW_KEY_4) renderer->incGPUTraversalMaxIters(10);
	if (key == GLFW_KEY_5) renderer->decPixelTolerance();
	if (key == GLFW_KEY_6) renderer->incPixelTolerance();
	if (key == GLFW_KEY_R) renderer->nextViewerRenderMode();
	if (key == GLFW_KEY_O) renderer->toggleUseMinDepthOptimization();

}

void updateSceneEvents(GLFWwindow *win) {

	// Custom scene updating events
	if (glfwGetKey(win, GLFW_KEY_DELETE) == GLFW_PRESS) {
		//if (_state == S_SVO || _state == S_DAG) {
		int nodeIndex = encoded_octree->getNodeIndex(cam->getCurrentConfig().pos);
		if (nodeIndex == -1) {
			nodeIndex = (std::rand() % encoded_octree->getNNodes());
		}
		renderer->clearVoxel(nodeIndex);
	}
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
	renderer->getCamera()->scroll_callback(xoffset, yoffset);
}

void updateInfo() {
	static char winMsg[256];
	const float drawTime = frameStats.time;
	const unsigned int raysSec = (unsigned int)(SCREEN_WIDTH * SCREEN_HEIGHT * 1000.0f  / drawTime);
	const float fps = 1000.0f / frameTime;
	sprintf(
		winMsg,
		"SymVox Viewer :: %.0f FPS | DrawT %.1f ms | %s Rays/s | Levls %i (%s^3) | MaxIters %i | PixTol %.2f | Beam Opt: %i  ::  http://vic.crs4.it",
		round(fps), drawTime,
		sl::human_readable_quantity(raysSec).c_str(),
		renderer->getDrawLevel(),
		sl::human_readable_quantity(1ULL << renderer->getDrawLevel(), "", 1024ULL).c_str(),
		renderer->getGPUTraversalMaxIters(),
		renderer->getPixelTolerance(),
		renderer->getUseMinDepthOptimization()
		);
	glfwSetWindowTitle(window, winMsg);
}

int main(int argc, char ** argv)
{
	
	printf(
		"\n"
		"====================================================================\n"
		"=====   SymVox Viewer (c)2016 ViC / CRS4   http://vic.crs4.it   ====\n"
		"====================================================================\n"
		"\n"
	);

	// Arguments parsing
	if (argc < 2) {
		printf("Usage: svviewer model.[svdag | ussvdag | ssvdag]");
		exit(1);
	}

	std::string inputFile(argv[1]);

	std::string ext = sl::pathname_extension(inputFile);
	bool incorrectFile = false;
	if (ext == "svdag" || ext == "SVDAG") {
		encoded_octree = new EncodedSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else if (ext == "ussvdag" || ext == "USSVDAG") {
		encoded_octree = new EncodedUSSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else if (ext == "ssvdag" || ext == "SSVDAG") {
		encoded_octree = new EncodedSSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else 
		incorrectFile = true;

	if (incorrectFile) {
		printf("* ERROR: Unsupported octree '%s'\n", inputFile.c_str());
		printf("         This viewer support '.svdag', '.ussvdag' and '.ssvdag' "
			"files, built with the SymVox tool 'svbuilder'.\n");
		return 1;
	}

	renderer = new OctreeDDARenderer(encoded_octree);
	renderer->setScreenResolution(SCREEN_WIDTH, SCREEN_HEIGHT);


	// Initialise GLFW --------------------------------
	if (!glfwInit()) {
		std::cerr << "ERROR: Failed to initialize GLFW" << std::endl;
		return -1;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Voxelator Viewer | Loading...", NULL, NULL);

	if (!window) {
		std::cerr << "ERROR: Can't open the window" << std::endl;
		glfwTerminate();
		return -1;
	}

	cam = new Camera(window);
	sl::aabox3f sceneBBox = renderer->getSceneBBox();
	cam->setInitCamera(sceneBBox.corner(6),
		sceneBBox.center(),
		Camera::Y_UP,
		60.0f,
		sceneBBox.diagonal().two_norm() * 0.001f,
		sceneBBox.diagonal().two_norm() * 10.0f);

	// Make the window's context current
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		printf("Failed to initialize GLEW\n");
		return -1;
	}

	// Hack to discard a bug in GLEW, that throuth an INVALID ENUM error
	// without much sense, when using core profile.
	// More info: https://www.opengl.org/wiki/OpenGL_Loading_Library
	glGetError();

	// Keyboard callback
	glfwSetKeyCallback(window, optionsKeyCallback);
	glfwSetScrollCallback(window, scrollCallback);

	renderer->setCamera(cam);
	renderer->init();
	renderer->selectRenderMode(OctreeDDARenderer::VIEWER);
	renderer->clearState();
	renderer->resetState();
	renderer->toggleRenderingStats();

	glfwSwapInterval(1); // Enable/disable vsync
	glClearColor(0.3f, 0.5, 0.7f, 1.0f);

	/////////// Render Loop //////////
	double t = glfwGetTime();
	unsigned int frameCount = 0, lastCountStep = 0;
	float et = 0;

	printHelp();

	do {
		// UPDATE ------------------
		cam->update();

		updateSceneEvents(window);

		// DRAW --------------------
		glClear(GL_COLOR_BUFFER_BIT);
		renderer->draw();
		frameCount++;

		glfwSwapBuffers(window);
		glfwPollEvents();

		// STATS ------------------
		et += float((glfwGetTime() - t)*1000.0);
		if (et > UPDATE_INFO_TIME) {
			frameStats = renderer->getStats();
			frameTime = et / float(frameCount - lastCountStep);
			updateInfo();
			lastCountStep = frameCount;
			et = 0;
		}
		t = glfwGetTime();

	} while (!glfwWindowShouldClose(window) && !finish);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
