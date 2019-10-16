
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

#include "octree_dda_renderer.hpp"
#include "halton.hpp"

#include <sl/clock.hpp>
#include <symvox/util.hpp>

#include <symvox/encoded_svdag.hpp>
#include <symvox/encoded_svo.hpp>
#include <symvox/encoded_ussvdag.hpp>
#include <symvox/encoded_ssvdag.hpp>

#define SHADER_FILE ("../shaders/octree_dda.frag.glsl")
#define NUM_MAX_HS_SAMPLES 256

OctreeDDARenderer::OctreeDDARenderer(EncodedOctree * eo) {

	_encodedOctree = eo;
	_selectedRenderMode = VIEWER;
	_shadowsInputDepthTex = 4;
	_shadowsInputPosTex = 5;
	_shadowsInputNormTex = 6;

	_selectedVoxelIndex = -1;
	_randomColors = true;
	_freqColors = true;
	_enableShadows = false;
	_lightPos = _encodedOctree->getSceneBBox()[1];

	if (dynamic_cast<EncodedSVO*>(_encodedOctree) != NULL) {
		_octreeFormat = SVO;
	}
	else if (dynamic_cast<EncodedSVDAG*>(_encodedOctree) != NULL) {
		_octreeFormat =SVDAG;
	}
	else if (dynamic_cast<EncodedUSSVDAG*>(_encodedOctree) != NULL) {
		_octreeFormat = USSVDAG;
	}
	else if (dynamic_cast<EncodedSSVDAG*>(_encodedOctree) != NULL) {
		_octreeFormat = SSVDAG;
	}
	else {
		printf("OctreeDDARenderer: ERROR! the EncodedOctree passed is not valid!\n");
		exit(-1);
	}

	_gpuTraversalMaxIters = 300;
	_drawLevel = _encodedOctree->getNLevels();
	_encodedOctree = eo;

	_fovy = 60.0f * 3.14159265359f / 180.0f;
	_pixelTolerance = 1.0;
	_drawLevel = 1;
	_viewerRenderMode = 0;

	_numAORays = 16;
	_lengthAORays = _encodedOctree->getHalfSide() / 100.0f;

	_sdag4UseTex3D = false;
}

void OctreeDDARenderer::uploadData(bool init) {
	GLint tmp;
	glGetIntegerv(GL_MAX_TEXTURE_BUFFER_SIZE, &tmp);
	const unsigned int maxTBOTexels = (unsigned int)tmp;
	glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, &tmp);
	const unsigned int maxT3DTexels = (unsigned int)tmp;
	const unsigned int maxT3DTexelsPow2 = maxT3DTexels * maxT3DTexels;

	// Clear data buffers in case they already were uploaded before
	// todo: ??
	//glDeleteTextures(4, _glTex);
	//glDeleteBuffers(2, _glBuf);

	if (_octreeFormat == SVO) {
		EncodedSVO * svo = static_cast<EncodedSVO *>(_encodedOctree);
		if (init) glGenBuffers(1, _glBuf);
		if (init) glGenTextures(1, _glTex);
		glBindTexture(GL_TEXTURE_BUFFER, _glTex[0]);
		glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[0]);
		if ((svo->getDataSize() / 16) > maxTBOTexels) printf("\t- WARNING! "); else printf("\t- ");
		printf("Uploading Data to TBO: %zu texels [%s]\n", svo->getDataSize() / 16, sl::human_readable_size(svo->getDataSize()).c_str());
		glBufferData(GL_TEXTURE_BUFFER, svo->getDataSize(), svo->getDataPtr(), GL_STATIC_DRAW);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, _glBuf[0]);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
		glBindBuffer(GL_TEXTURE_BUFFER, 0);
	}
	else if (_octreeFormat == SVDAG) {
		EncodedSVDAG * dag = static_cast<EncodedSVDAG *>(_encodedOctree);
		if (init) glGenBuffers(1, _glBuf);
		if (init) glGenTextures(1, _glTex);
		glBindTexture(GL_TEXTURE_BUFFER, _glTex[0]);
		glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[0]);
		if ((dag->getDataSize() / 16) > maxTBOTexels) printf("\t- WARNING! "); else printf("\t- ");
		printf("Uploading Data to TBO: %zu texels [%s]\n", dag->getDataSize() / 16, sl::human_readable_size(dag->getDataSize()).c_str());
		glBufferData(GL_TEXTURE_BUFFER, dag->getDataSize(), dag->getDataPtr(), GL_DYNAMIC_DRAW);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, _glBuf[0]);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
		glBindBuffer(GL_TEXTURE_BUFFER, 0);
	}
	else if (_octreeFormat == USSVDAG) {
		EncodedUSSVDAG * sdag1 = static_cast<EncodedUSSVDAG *>(_encodedOctree);
		if (init) glGenBuffers(1, _glBuf);
		if (init) glGenTextures(1, _glTex);
		glBindTexture(GL_TEXTURE_BUFFER, _glTex[0]);
		glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[0]);
		if ((sdag1->getDataSize() / 16) > maxTBOTexels) printf("\t- WARNING! "); else printf("\t- ");
		printf("Uploading Data to TBO: %zu texels [%s]\n", sdag1->getDataSize() / 16, sl::human_readable_size(sdag1->getDataSize()).c_str());
		glBufferData(GL_TEXTURE_BUFFER, sdag1->getDataSize(), sdag1->getDataPtr(), GL_STATIC_DRAW);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, _glBuf[0]);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
		glBindBuffer(GL_TEXTURE_BUFFER, 0);
	}
	else if (_octreeFormat == SSVDAG) {
		EncodedSSVDAG * sdag4 = static_cast<EncodedSSVDAG *>(_encodedOctree);

		// TEXTURE_0 -> TBO with LEAVES
		if (init) glGenBuffers(1, _glBuf);
		if (init) glGenTextures(1, &_glTex[0]);
		glBindTexture(GL_TEXTURE_BUFFER, _glTex[0]);
		glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[0]);
		if ((sdag4->getLeafNodesSize() / 8) > maxTBOTexels) printf("\t- WARNING! "); else printf("\t- ");
		printf("Uploading Leaf nodes to TBO: %u texels [%s]\n", sdag4->getLeafNodesSize() / 8, sl::human_readable_size(sdag4->getLeafNodesSize()).c_str());
		glBufferData(GL_TEXTURE_BUFFER, sdag4->getLeafNodesSize(), sdag4->getLeafNodesDataPtr(), GL_STATIC_DRAW);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32UI, _glBuf[0]);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
		glBindBuffer(GL_TEXTURE_BUFFER, 0);

		// TEXTURE_1 -> Inners
		if (init) glGenTextures(1, &_glTex[1]);
		if (_sdag4UseTex3D) { // as a 3D Texture
			const int texelSize = 2; // corresponding to GL_R16UI
			const unsigned int neededTexels = sdag4->getNInnerNodes();
			unsigned int depthLayers = 1 + (neededTexels / maxT3DTexelsPow2);
			const unsigned int texelsToAllocate = maxT3DTexelsPow2 * depthLayers;
			sdag4->expandInnerBuffer(texelsToAllocate);
			glBindTexture(GL_TEXTURE_3D, _glTex[1]);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
			printf("\t -Uploading Inners to 3D texture (%i x %i x %i) [%s]\n", maxT3DTexels, maxT3DTexels, depthLayers, sl::human_readable_size(maxT3DTexelsPow2 * depthLayers * 2).c_str());
			glTexImage3D(GL_TEXTURE_3D, 0, GL_R16UI, maxT3DTexels, maxT3DTexels, depthLayers, 0, GL_RED_INTEGER, GL_UNSIGNED_SHORT, sdag4->getInnerNodesDataPtr());
			glBindTexture(GL_TEXTURE_3D, 0);
			sdag4->compressInnerBuffer(neededTexels);
		}
		else { // as a TBO
			if (init) glGenBuffers(1, &_glBuf[1]);
			glBindTexture(GL_TEXTURE_BUFFER, _glTex[1]);
			glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[1]);
			if ((sdag4->getInnerNodesSize() / 8) > maxTBOTexels) printf("\t- WARNING! "); else printf("\t- ");
			printf("Uploading Inner nodes to TBO: %u texels [%s]\n", sdag4->getInnerNodesSize() / 8, sl::human_readable_size(sdag4->getInnerNodesSize()).c_str());
			glBufferData(GL_TEXTURE_BUFFER, sdag4->getInnerNodesSize(), sdag4->getInnerNodesDataPtr(), GL_STATIC_DRAW);
			glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA16UI, _glBuf[1]);
			glBindTexture(GL_TEXTURE_BUFFER, 0);
			glBindBuffer(GL_TEXTURE_BUFFER, 0);
		}
		glUniform1uiv(
			_program[_selectedRenderMode].getUniformID("levelOffsets"),
			sdag4->getNLevelOffsets(),
			(GLuint *)sdag4->getLevelOffsetsPtr());
		
		glUseProgram(_program[DEPTH].getGLProgram());
		glUniform1uiv(
			_program[DEPTH].getUniformID("levelOffsets"),
			sdag4->getNLevelOffsets(),
			(GLuint *)sdag4->getLevelOffsetsPtr());
		
	}

	if (!init) {
		for (int i = 0; i < 4; ++i) {
			// Todo: Properly reset uniforms when reuploading data
			// Just copied some lines for now, mostly for SVDAG
			_program[i].setDefine("LEVELS", (int)_encodedOctree->getNLevels());
			_program[i].setDefine("ROOT_HALF_SIDE", _encodedOctree->getRootSide() / 2.0f);

			if (_octreeFormat == SSVDAG) {
				_program[i].setDefine("INNER_LEVELS", (int)(_encodedOctree->getNLevels() - 2));
			}
			else {
				_program[i].setDefine("INNER_LEVELS", (int)(_encodedOctree->getNLevels() - 1));
			}

			// Non-dynamic generic uniforms
			const sl::aabox3f& box = _encodedOctree->getSceneBBox();
			glUseProgram(_program[i].getGLProgram());
			glUniform3fv(_program[i].getUniformID("sceneBBoxMin"), 1, box[0].to_pointer());
			glUniform3fv(_program[i].getUniformID("sceneBBoxMax"), 1, box[1].to_pointer());
			glUniform3fv(_program[i].getUniformID("sceneCenter"), 1, box.center().to_pointer());
			glUniform1f(_program[i].getUniformID("rootHalfSide"), _encodedOctree->getHalfSide());
		}
	}
}

void OctreeDDARenderer::init() {

	_glTex[2] = generateChildrenIndirectionTex();

	printf("* Initializing Octree DDA Renderer in a %s\n", getOctreeFormatName().c_str());

	_drawLevel = _encodedOctree->getNLevels();
	_projectionFactor = getProjectionFactor(_pixelTolerance);
	generateHSSamples(NUM_MAX_HS_SAMPLES);

	GLint tmp;
	glGetIntegerv(GL_MAX_TEXTURE_BUFFER_SIZE, &tmp);
	const unsigned int maxTBOTexels = (unsigned int)tmp;
	glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, &tmp);
	const unsigned int maxT3DTexels = (unsigned int)tmp;
	const unsigned int maxT3DTexelsPow2 = maxT3DTexels * maxT3DTexels;

	printf("\t- LIMITS: Max TBO %i texels - Max Tex3D %i texels\n", maxTBOTexels, maxT3DTexels);

	if (_octreeFormat == SSVDAG) {
		sl::uint32_t maxInnerIndirection = 4 * maxTBOTexels * 2; //  RGBA16UI
		if (static_cast<EncodedSSVDAG *>(_encodedOctree)->getInnerNodesSize() > maxInnerIndirection)
			_sdag4UseTex3D = true;
	}

	// create the GLSL programs (VIEWER, DEPTH, SHADOW)
	for (int i = 0; i < 4; ++i) {
		std::string name = "Octree DDA Renderer [" + getOctreeFormatName() + ", " + getRenderModeName(RenderMode(i)) + "]";

		// Setup GLSL Program //////////////////////////////////////////////////////
		printf("\t- Setup and compile/link [%s, %s] GLSL program... \n",
			getRenderModeName(RenderMode(i)).c_str(),
			getOctreeFormatName().c_str());

		_program[i].setName(name);

		if (RenderMode(i) == VIEWER) _program[i].setDefine("VIEWER_MODE", 1);
		else if (RenderMode(i) == DEPTH) _program[i].setDefine("DEPTH_MODE", 1);
		else if (RenderMode(i) == SHADOW) _program[i].setDefine("SHADOW_MODE", 1);
		else if (RenderMode(i) == AO) {
			_program[i].setDefine("AO_MODE", 1);
			_program[i].setDefine("N_HS_SAMPLES", NUM_MAX_HS_SAMPLES);
		}

		_program[i].setDefine("LEVELS", (int)_encodedOctree->getNLevels());
		_program[i].setDefine("INNER_LEVELS", (int)(_encodedOctree->getNLevels() - 2));
		_program[i].setDefine("ROOT_HALF_SIDE", _encodedOctree->getRootSide() / 2.0f);
		if (_octreeFormat == SSVDAG) {
			_program[i].setDefine("INNER_LEVELS", (int)(_encodedOctree->getNLevels() - 2));
			if (_sdag4UseTex3D) {
				_program[i].setDefine("SSVDAG_TEX3D", 1);
				_program[i].setDefine("TEX3D_SIZE", maxT3DTexels);
				_program[i].setDefine("TEX3D_SIZE_POW2", maxT3DTexelsPow2);
			}
		}
		else {
			_program[i].setDefine("INNER_LEVELS", (int)(_encodedOctree->getNLevels() - 1));
		}
		switch (_octreeFormat) {
		case SVO:
			_program[i].setDefine("SVO", 1); break;
		case SVDAG:
			_program[i].setDefine("SVDAG", 1); break;
		case USSVDAG:
			_program[i].setDefine("USSVDAG", 1); break;
		case SSVDAG:
			_program[i].setDefine("SSVDAG", 1); break;
		default:
			return;
		}

		printf("\t\t - Added [%i] lines to the start of the shader program\n", _program[i].getNumDefines() + 1);

		if (!_program[i].setAndCompileShader(GLSLProgram::VERTEX_SHADER, _sqr.getDefaultVertexProgram())) exit(-1);
		if (!_program[i].loadAndCompileShader(GLSLProgram::FRAGMENT_SHADER, SHADER_FILE)) exit(-1);
		if (!_program[i].link()) exit(-1);

		// Non-dynamic generic uniforms
		const sl::aabox3f& box = _encodedOctree->getSceneBBox();
		glUseProgram(_program[i].getGLProgram());
		glUniform3fv(_program[i].getUniformID("sceneBBoxMin"), 1, box[0].to_pointer());
		glUniform3fv(_program[i].getUniformID("sceneBBoxMax"), 1, box[1].to_pointer());
		glUniform3fv(_program[i].getUniformID("sceneCenter"), 1, box.center().to_pointer());
		glUniform1f(_program[i].getUniformID("rootHalfSide"), _encodedOctree->getHalfSide());
		if (RenderMode(i) == SHADOW || RenderMode(i) == AO) {
			glUniform1i(_program[i].getUniformID("depthTex"), 4);
			glUniform1i(_program[i].getUniformID("hitPosTex"), 5);
			glUniform1i(_program[i].getUniformID("hitNormTex"), 6);
		}
		else {
			glUniform2fv(_program[i].getUniformID("screenRes"), 1, _screenRes.to_pointer());
			glUniform1i(_program[i].getUniformID("minDepthTex"), 3);
		}
		if (RenderMode(i) == AO) {
			glUniform3fv(_program[i].getUniformID("hsSamples"), (GLsizei)_hsSamples3D.size(), (GLfloat *)&_hsSamples3D[0]);
		}

		if (_octreeFormat == SSVDAG) {
			EncodedSSVDAG * sdag4 = static_cast<EncodedSSVDAG *>(_encodedOctree);
			glUniform1i(_program[i].getUniformID("leafNodes"), 0);
			glUniform1i(_program[i].getUniformID("innerNodes"), 1);
			glUniform1uiv(_program[i].getUniformID("levelOffsets"), sdag4->getNLevelOffsets(), (GLuint *)sdag4->getLevelOffsetsPtr());
			glUniform1i(_program[i].getUniformID("childIndir"), 2);
		} else {
			glUniform1i(_program[i].getUniformID("nodes"), 0);
		}

		printf("\t\t - OK!\n");
		printfGLError();
	}

	uploadData(true);

	glGenFramebuffers(1, &_glFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, _glFBO);
	glGenTextures(1, &_glTex[3]);
	glBindTexture(GL_TEXTURE_2D, _glTex[3]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, (int)_screenRes[0]/8, (int)_screenRes[1]/8, 0, GL_RGB, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, _glTex[3], 0);
	GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, DrawBuffers);
	glBindTexture(GL_TEXTURE_2D, 0);

	checkCurrentFrameBuffer();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	printfGLError();

	_sqr.initGL();

	resetState();
}

void OctreeDDARenderer::clearState() {
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_BUFFER, 0);
	if (_octreeFormat == SSVDAG) {
		glActiveTexture(GL_TEXTURE1);
		glBindTexture((_sdag4UseTex3D ? GL_TEXTURE_3D : GL_TEXTURE_BUFFER), 0);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER, 0);
	}
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void OctreeDDARenderer::resetState() {
	glDisable(GL_DEPTH_TEST);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_BUFFER, _glTex[0]);

	if (_octreeFormat == SSVDAG) {
		glActiveTexture(GL_TEXTURE1);
		glBindTexture((_sdag4UseTex3D ? GL_TEXTURE_3D : GL_TEXTURE_BUFFER), _glTex[1]);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER, _glTex[2]);
	}

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, _glTex[3]);
}

void OctreeDDARenderer::draw(GLint renderBuffer) {

	//resetState();
	_rendererMonitor.beginFrame();
	printfGLError();

	if (_useMinDepthOptimization && (_selectedRenderMode == VIEWER || _selectedRenderMode == DEPTH)) {
		glUseProgram(_program[DEPTH].getGLProgram());
		glUniform1ui(_program[DEPTH].getUniformID("maxIters"), _gpuTraversalMaxIters);
		glUniform1ui(_program[DEPTH].getUniformID("drawLevel"), _drawLevel);
		glUniform1i(_program[DEPTH].getUniformID("useMinDepthTex"), false);
		glUniformMatrix4fv(_program[DEPTH].getUniformID("viewMatInv"), 1, GL_FALSE, _camera->getViewMatrixInv().to_pointer());
		glUniformMatrix4fv(_program[DEPTH].getUniformID("projMatInv"), 1, GL_FALSE, _camera->getProjMatrixInv().to_pointer());
		_projectionFactor = getProjectionFactor(1.0f, 8.0f); // FIXME ? ? ? ? ? ? ? ? ? ? ? ? ? ? ?
		glUniform1f(_program[DEPTH].getUniformID("projectionFactor"), _projectionFactor);
		glUniform2f(_program[DEPTH].getUniformID("screenRes"), _screenRes[0] / 8, _screenRes[1] / 8);
		glBindFramebuffer(GL_FRAMEBUFFER, _glFBO);
		glClearColor(0.0, 0.0, 0.0, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		_sqr.draw();
		glBindFramebuffer(GL_FRAMEBUFFER, renderBuffer);
		glUniform2f(_program[DEPTH].getUniformID("screenRes"), _screenRes[0], _screenRes[1]);
		_projectionFactor = getProjectionFactor(_pixelTolerance);
		glUseProgram(_program[_selectedRenderMode].getGLProgram());
		glUniform1i(_program[_selectedRenderMode].getUniformID("useMinDepthTex"), true);
	}
	else {
		glUseProgram(_program[_selectedRenderMode].getGLProgram());
		if (_selectedRenderMode == VIEWER || _selectedRenderMode == DEPTH) glUniform1i(_program[_selectedRenderMode].getUniformID("useMinDepthTex"), false);
	}

	//resetState();
	glUseProgram(_program[_selectedRenderMode].getGLProgram());
	glUniform1ui(_program[_selectedRenderMode].getUniformID("maxIters"), _gpuTraversalMaxIters);
	glUniform1ui(_program[_selectedRenderMode].getUniformID("drawLevel"), _drawLevel);

	if (_selectedRenderMode != SHADOW && _selectedRenderMode != AO) {
		glUniformMatrix4fv(_program[_selectedRenderMode].getUniformID("viewMatInv"), 1, GL_FALSE, _camera->getViewMatrixInv().to_pointer());
		glUniformMatrix4fv(_program[_selectedRenderMode].getUniformID("projMatInv"), 1, GL_FALSE, _camera->getProjMatrixInv().to_pointer());
		glUniform1f(_program[_selectedRenderMode].getUniformID("projectionFactor"), _projectionFactor);
	}

	if (_selectedRenderMode == VIEWER) {
		glUniform1i(_program[_selectedRenderMode].getUniformID("viewerRenderMode"), _viewerRenderMode);
		glUniform1ui(_program[_selectedRenderMode].getUniformID("selectedVoxelIndex"), _selectedVoxelIndex);
		glUniform1i(_program[_selectedRenderMode].getUniformID("randomColors"), _randomColors);
		glUniform1i(_program[_selectedRenderMode].getUniformID("enableShadows"), _enableShadows);
		if (_octreeFormat == SSVDAG) glUniform1i(_program[_selectedRenderMode].getUniformID("freqColors"), _freqColors);
	}
	if (_selectedRenderMode == SHADOW || _selectedRenderMode == VIEWER)
		glUniform3fv(_program[_selectedRenderMode].getUniformID("lightPos"), 1, _lightPos.to_pointer());

	if (_selectedRenderMode == AO) {
		glUniform1i(_program[_selectedRenderMode].getUniformID("numAORays"), _numAORays);
		glUniform1f(_program[_selectedRenderMode].getUniformID("lengthAORays"), _lengthAORays);
	}

	_sqr.draw();
	_rendererMonitor.endFrame();
	printfGLError();
	//clearState();
}

void OctreeDDARenderer::end() {
	clearState();
	glDeleteTextures(4, _glTex);
	glDeleteBuffers(2, _glBuf);
}

void OctreeDDARenderer::clearVoxel(int position)
{
	EncodedSVDAG * dag = static_cast<EncodedSVDAG *>(_encodedOctree);

	int numNodes = (dag->getDataSize() / sizeof(sl::uint32_t));
	int nodeIndex = position;

	if (nodeIndex < 0 || nodeIndex > numNodes) {
		printf("\t- Voxel index for deletion is out of bounds! %u / %u \n", nodeIndex, numNodes);
		return;
	}
	else {
		printf("\t- Deleting voxel index %u \n", nodeIndex);
	}

	printf("Index: %u, sizeof uint32_t: %zu, uint32_t(0): %u \n", nodeIndex, sizeof(sl::uint32_t), sl::uint32_t(0));

    // Instead of replacing with 0,
    // copy leaf to end of node list, update pointer to it
	sl::uint32_t subData[] = { sl::uint32_t(0) };

	glBindBuffer(GL_TEXTURE_BUFFER, _glBuf[0]);
	glBufferSubData(GL_TEXTURE_BUFFER, nodeIndex * sizeof(sl::uint32_t), sizeof(sl::uint32_t), &subData);
	glBindBuffer(GL_TEXTURE_BUFFER, 0);
}

std::string OctreeDDARenderer::getOctreeFormatName() {
	std::string s;
	switch (_octreeFormat) {
	case SVO: s = "SVO"; break;
	case SVDAG: s = "SVDAG"; break;
	case USSVDAG: s = "USSVDAG"; break;
	case SSVDAG: s = "SSVDAG"; break;
	default: s = "NO_VALID_OCTREE_FORMAT";
	}
	return s;
}

std::string OctreeDDARenderer::getRenderModeName(RenderMode mode) {
	std::string s;
	switch (mode) {
	case VIEWER: s = "VIEWER"; break;
	case DEPTH: s = "DEPTH"; break;
	case SHADOW: s = "SHADOW"; break;
	case AO: s = "AO"; break;
	default: s = "NO_VALID_RENDER_MODE";
	}
	return s;
}

float OctreeDDARenderer::getProjectionFactor(float pixelTolerance, float screenDivisor) {
	const float inv_2tan_half_fovy = 1.0f / (2.0f * tan(0.5f * _fovy));
	const float screen_tolerance = _pixelTolerance / (_screenRes[1]/ screenDivisor);
	return inv_2tan_half_fovy / screen_tolerance;
}


void OctreeDDARenderer::generateHSSamples(int n) {
	_hsSamples2D.reserve(n);
	_hsSamples3D.reserve(n);

	HaltonSequence hs(2);

	for (int i = 0; i < n; i++) {
		sl::vector2f s2;
		hs.getNext(s2.to_pointer());

		_hsSamples2D.push_back(s2 - sl::vector2f(0.5f, 0.5f));

		// for hemisphere: distribute samples according to cosine weighted solid angle
		sl::vector3f s3;
		s3[1] = sqrt(s2[0]);
		const float l = sqrt(1.0f - s3[1] * s3[1]);
		s3[0] = l * cos(2.0f * s2[1] * 3.141592654f);
		s3[2] = l * sin(2.0f * s2[1] * 3.141592654f);
		_hsSamples3D.push_back(s3);
		//printf("%.2f\t%.2f\t%.2f\n", s3[0], s3[1], s3[2]);
	}
}

GLuint OctreeDDARenderer::generateChildrenIndirectionTex() {

	sl::int8_t data[65536 * 8];

	for (sl::uint32_t i = 0; i < 65536; ++i) {
		//printf("%i\t", i);
		for (unsigned int j = 0; j < 8; ++j) {
			sl::int8_t val = 0;
			for (unsigned int k = 7; k > j; --k) {
				val += sl::min(sl::int8_t((i >> (2 * k)) & 3), sl::int8_t(2));
			}
			data[j * 65536 + i] = val;
			//printf("%i ", val);
		}
		//printf("\n");
	}

	GLuint tex, buf;
	glGenBuffers(1, &buf);
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_BUFFER, tex);
	glBindBuffer(GL_TEXTURE_BUFFER, buf);
	glBufferData(GL_TEXTURE_BUFFER, 65536 * 8, (GLvoid *)data, GL_STATIC_DRAW);
	glTexBuffer(GL_TEXTURE_BUFFER, GL_R8I, buf);
	glBindTexture(GL_TEXTURE_BUFFER, 0);
	glBindBuffer(GL_TEXTURE_BUFFER, 0);

	return tex;
}
