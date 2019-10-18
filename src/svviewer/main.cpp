
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

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#define UPDATE_INFO_TIME 200.0 // ms
#define SCREEN_WIDTH	1280
#define SCREEN_HEIGHT	720

bool finish = false;
bool printCamera = false;
float frameTime = 0;
std::string filename = "";
static char filenameInput[128] = "";
static char filenameInput2[128] = "";

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
    ImGuiIO& io = ImGui::GetIO();
	if (action != GLFW_PRESS || io.WantCaptureKeyboard) return;

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
	if (key == GLFW_KEY_K && mods == 0) renderer->toggleRandomColors();
	if (key == GLFW_KEY_K && mods == GLFW_MOD_CONTROL) renderer->toggleFreqColors();
	if (key == GLFW_KEY_L) renderer->setLightPos(renderer->getCamera()->getCurrentConfig().pos);
	if (key == GLFW_KEY_F) cam->toggleCamController();

	if (key == GLFW_KEY_F1) renderer->selectRenderMode(OctreeDDARenderer::RenderMode::AO);
	if (key == GLFW_KEY_F2) renderer->selectRenderMode(OctreeDDARenderer::RenderMode::DEPTH);
	if (key == GLFW_KEY_F3) renderer->selectRenderMode(OctreeDDARenderer::RenderMode::SHADOW);
	if (key == GLFW_KEY_F4) renderer->selectRenderMode(OctreeDDARenderer::RenderMode::VIEWER);
}

static sl::vector3f fromHomog(const sl::vector4f v) { return sl::vector3f(v[0] / v[3], v[1] / v[3], v[2] / v[3]); }
static sl::vector3f v4ToV3(const sl::vector4f v) { return sl::vector3f(v[0], v[1], v[2]); }

/** Returns -1 when no voxel is found */
int getVoxelIndexAtCursor() {
	// Transform [0, 1] cursor coords to [-1, 1]
	sl::point2f coords = cam->getCursorCoordsNorm();
	coords[0] = coords[0] * 2 - 1;
	coords[1] = 1.f - coords[1] * 2;
	printf("mouse coords: %f, %f - \t", coords[0], coords[1]);

	// Apply projection matrix
	sl::vector3f rayStart = fromHomog(cam->getProjMatrixInv().as_matrix() * sl::vector4f(coords[0], coords[1], 0, 1));
	sl::vector3f rayEnd = fromHomog(cam->getProjMatrixInv().as_matrix() * sl::vector4f(coords[0], coords[1], 1, 1));
	sl::vector3f rayDir = (rayEnd - rayStart).ok_normalized();
	// Apply view matrix
	rayStart = v4ToV3(cam->getViewMatrixInv().as_matrix() * sl::vector4f(rayStart[0], rayStart[1], rayStart[2], 1));
	rayEnd = v4ToV3(cam->getViewMatrixInv().as_matrix() * sl::vector4f(rayEnd[0], rayEnd[1], rayEnd[2], 1));
	rayDir = (rayEnd - rayStart).ok_normalized();
	printf("mouse ray dir: %f, %f, %f \n", rayDir[0], rayDir[1], rayDir[2]);

	sl::point3f delPos = cam->getCurrentConfig().pos;
	for (int i = 0; i < 1024; i++) {
		// Traverse ray along small increments.
		// Todo: use stack + dda
		float stepSize = encoded_octree->getHalfSide(renderer->getDrawLevel());
		delPos += rayDir * stepSize;

		int nodeIndex = encoded_octree->getNodeIndex(delPos, renderer->getDrawLevel());
		if (nodeIndex == -1) {
			// renderer->setSelectedVoxelIndex(0);
			continue;
			//return;
			// nodeIndex = (std::rand() % encoded_octree->getNNodes());
		}
		else {
			// printf("\nDELETED A VOXEL!!!!!!!!!!!!!!!!!!!! \n");
			// break;
			return nodeIndex;
		}
	}
	return -1;
}

// Custom scene updating events
void updateSceneEvents(GLFWwindow *win) {

	bool doHighlight = glfwGetKey(win, GLFW_KEY_ENTER) == GLFW_PRESS;
	bool doDelete = glfwGetKey(win, GLFW_KEY_DELETE) == GLFW_PRESS;
	
	if (doHighlight || doDelete) {
		//if (_state == S_SVO || _state == S_DAG) {

		int nodeIndex = getVoxelIndexAtCursor();

		if (doHighlight) {
			renderer->setSelectedVoxelIndex(nodeIndex == -1 ? 0 : nodeIndex);
		}
		else if (doDelete) {
			renderer->clearVoxel(nodeIndex);
		}
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

bool loadFile(std::string inputFile) {
	filename = sl::pathname_base(inputFile);
	std::string ext = sl::pathname_extension(inputFile);
	bool incorrectFile = false;

	if (encoded_octree != nullptr)
        delete encoded_octree;

	if (ext == "svdag" || ext == "SVDAG") {
		encoded_octree = new EncodedSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else if (ext == "ussvdag" || ext == "USSVDAG") {
		encoded_octree = new EncodedUSSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else if (ext == "ssvdag" || ext == "SSVDAG"
		  || ext == "esvdag" || ext == "ESVDAG") {
		encoded_octree = new EncodedSSVDAG();
		if (!encoded_octree->load(inputFile)) incorrectFile = true;
	}
	else
		incorrectFile = true;
	return incorrectFile;
}

void handleImgui() {

	int drawLevelInput = renderer->getDrawLevel();
	int maxItersInput = renderer->getGPUTraversalMaxIters();
	bool beamOptInput = renderer->getUseMinDepthOptimization();
	bool randomFreqInput = renderer->getFreqColors();
	bool randomColsInput = renderer->getRandomColors();
	bool shadowsEnabledInput = renderer->getShadowsEnabled();
	int normSamplesInput = renderer->getNormSamples();
	bool isUsingOrbitController = renderer->getCamera()->isUsingOrbitController();
	float pixTolInput = renderer->getPixelTolerance();
	float fovInput = renderer->getFovH();
	float walkFactorInput = renderer->getCamera()->getWalkFactor();

	const sl::point3f &camPos = renderer->getCamera()->getCurrentConfig().pos;
	float camPosInput[]{ 0,0,0 };
	camPosInput[0] = camPos[0]; camPosInput[1] = camPos[1]; camPosInput[2] = camPos[2];

	const sl::point3f &lightPos = renderer->getLightPos();
	float lightPosInput[]{ 0,0,0 };
	lightPosInput[0] = lightPos[0]; lightPosInput[1] = lightPos[1]; lightPosInput[2] = lightPos[2];

	const static char* renderModes[] = { "ITERATIONS", "DEPTH", "LEVELS", "PRETTY" };
	int renderModeInput = renderer->getViewerRenderMode();

	if (ImGui::Begin("SymVox - Fork by RvanderLaan")) {

        ImGui::Text("File: \"%s\" - %s", filename.c_str(), encoded_octree->getDescription().c_str());

        bool doLoadFile1 = ImGui::InputText("##file-slot-1", filenameInput, IM_ARRAYSIZE(filenameInput),
                                            ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::SameLine();
        doLoadFile1 = doLoadFile1 || ImGui::Button("Load file slot 1");
        if (doLoadFile1) {
            bool error = loadFile(filenameInput);
            if (!error) {
                renderer->setEncodedOctree(encoded_octree);
                renderer->uploadData();
            }
        }
        bool doLoadFile2 = ImGui::InputText("##file-slot-2", filenameInput2, IM_ARRAYSIZE(filenameInput2),
                                            ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::SameLine();
        doLoadFile2 = doLoadFile2 || ImGui::Button("Load file slot 2");
        if (doLoadFile2) {
            bool error = loadFile(filenameInput2);
            if (!error) {
                renderer->setEncodedOctree(encoded_octree);
                renderer->uploadData();
            }
        }

        if (ImGui::SliderInt("Draw level", &drawLevelInput, 1, encoded_octree->getNLevels()))
            renderer->setDrawLevel(drawLevelInput);
        if (ImGui::SliderInt("Max trav iters", &maxItersInput, 1, 500))
            renderer->setGPUTraversalMaxIters(maxItersInput);
        if (ImGui::SliderFloat("Pixel tolerance", &pixTolInput, 0.0f, 4.0f))
            renderer->setPixelTolerance(pixTolInput);
        if (ImGui::Checkbox("Beam optimization", &beamOptInput))
            renderer->toggleUseMinDepthOptimization();
		if (ImGui::Checkbox("Frequency colors (E/SSVDAG)", &randomFreqInput))
            renderer->toggleFreqColors();
        if (ImGui::Checkbox("Random colors", &randomColsInput))
            renderer->toggleRandomColors();
		if (ImGui::Checkbox("Shadows enabled ", &shadowsEnabledInput))
			renderer->toggleShadowsEnabled();
		if (ImGui::SliderInt("Normal samples ", &normSamplesInput, 0, 16))
            renderer->setNormSamples(normSamplesInput);
       
        //if (ImGui::SliderFloat("Fov", &fovInput, 0.0f, 6.28))
        //	renderer->setFovH(fovInput);


		ImGui::Text("Camera controller");
		ImGui::SameLine();
        if (ImGui::RadioButton("Orbit", isUsingOrbitController))
            renderer->getCamera()->toggleCamController();
		ImGui::SameLine();
		if (ImGui::RadioButton("First person", !isUsingOrbitController))
            renderer->getCamera()->toggleCamController();

        if (ImGui::SliderFloat("Move speed", &walkFactorInput, 0.0f, 4.0f))
            renderer->getCamera()->setWalkFactor(walkFactorInput);

        if (ImGui::DragFloat3("Camera position", camPosInput, 2)) {
			sl::vector3f delta(camPosInput[0], camPosInput[1], camPosInput[2]);
			delta -= renderer->getCamera()->getCurrentConfig().pos.as_vector();
			renderer->getCamera()->getCurrentConfig().pos += delta;
			renderer->getCamera()->getCurrentConfig().target += delta;
			renderer->getCamera()->updateMatrices();
        }

		if (ImGui::DragFloat3("Light position", lightPosInput, 2)) {
			renderer->setLightPos(sl::point3f(lightPosInput[0], lightPosInput[1], lightPosInput[2]));
		}

        if (ImGui::Combo("Render mode", &renderModeInput, renderModes, IM_ARRAYSIZE(renderModes))) {
            renderer->setViewerRenderMode(renderModeInput);
        }

        if (ImGui::Button("Set light pos at camera")) {
            renderer->setLightPos(renderer->getCamera()->getCurrentConfig().pos);
        }

        // todo: show count, # references (per level), etc. + delete button
        ImGui::Text("Selected voxel index (Hover + Enter): %i", renderer->getSelectedVoxelIndex());
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate,
                    ImGui::GetIO().Framerate);
    }
	ImGui::End();
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

	// Initialize imgui text inputs
    strcpy(filenameInput, inputFile.c_str());
    strcpy(filenameInput2, inputFile.c_str());

	std::string ext = sl::pathname_extension(inputFile);
	bool incorrectFile = loadFile(inputFile);

	if (incorrectFile) {
		printf("* ERROR: Unsupported octree '%s'\n", inputFile.c_str());
		printf("         This viewer support '.svdag', '.ussvdag' and '.ssvdag' "
			"files, built with the SymVox tool 'svbuilder'.\n");
		return 1;
	}

//    encoded_octree->getRootTravNode()
//    return 0;

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
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

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

	cam->setWalkFactor(sceneBBox.diagonal().two_norm() * 0.001f);

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

	/////////// imgui setup //////////
	const char* glsl_version = "#version 440";
	ImGui::CreateContext();
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	/////////// Render Loop //////////
	double t = glfwGetTime();
	unsigned int frameCount = 0, lastCountStep = 0;
	float et = 0;

	printHelp();

	do {
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// UPDATE ------------------
		cam->update();

		updateSceneEvents(window);

		// DRAW --------------------
		glClear(GL_COLOR_BUFFER_BIT);
		renderer->draw();
		frameCount++;

		handleImgui();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

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

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}
